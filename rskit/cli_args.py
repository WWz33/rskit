import shlex
from typing import Iterable, List, Sequence, Set


def split_extra_args(extra_args: str) -> List[str]:
    """Split a user-provided passthrough argument string."""
    if not extra_args:
        return []

    lexer = shlex.shlex(extra_args, posix=True)
    lexer.whitespace_split = True
    lexer.escape = ""
    return list(lexer)


def merge_extra_args(
    default_args: Sequence[str],
    extra_args: str,
    protected_options: Iterable[str],
) -> List[str]:
    """Merge passthrough CLI args, letting user-provided options replace defaults."""
    user_args = split_extra_args(extra_args)
    if not user_args:
        return list(default_args)

    protected = set(protected_options)
    user_options = _canonical_options(_option_names(user_args), protected)
    blocked = sorted(user_options & protected)
    if blocked:
        blocked_text = ", ".join(blocked)
        raise ValueError(f"Protected options cannot be set through passthrough args: {blocked_text}")

    merged = _remove_replaced_defaults(list(default_args), user_options, protected)
    merged.extend(user_args)
    return merged


def _option_names(args: Sequence[str]) -> Set[str]:
    return {_option_name(arg) for arg in args if _is_option(arg)}


def _canonical_options(options: Set[str], protected: Set[str]) -> Set[str]:
    """Expand attached short values and long-option prefixes to protected full names.

    Target tools accept them (salmon/boost guess unambiguous prefixes, short flags
    accept attached values like `-p8`), so the guard must recognize them too.
    """
    canonical = set(options)
    for option in options:
        if option in protected:
            continue
        if not option.startswith("--"):
            token = option[:2] if len(option) > 2 else option
            if token in protected:
                canonical.add(token)
            continue
        canonical.update(p for p in protected if p.startswith(option))
    return canonical


def _remove_replaced_defaults(
    default_args: List[str],
    user_options: Set[str],
    protected_options: Set[str],
) -> List[str]:
    merged = []
    index = 0
    while index < len(default_args):
        arg = default_args[index]
        if arg in user_options and arg not in protected_options:
            index += 1
            while index < len(default_args) and not _is_option(default_args[index]):
                index += 1
            continue

        merged.append(arg)
        index += 1

    return merged


def _is_option(arg: str) -> bool:
    return arg.startswith("-") and arg != "-"


def _option_name(arg: str) -> str:
    return arg.split("=", 1)[0]
