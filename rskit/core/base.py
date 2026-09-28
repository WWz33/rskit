import shutil
import subprocess
from abc import ABC, abstractmethod
from typing import Optional
from rskit.utils.logger import get_logger


def require_tools(*tool_names: str) -> None:
    """Fail fast when required external tools are missing from PATH.

    Called before any expensive work so the user sees a clear message instead
    of a FileNotFoundError traceback from deep inside a subprocess call.
    """
    missing = [name for name in tool_names if shutil.which(name) is None]
    if missing:
        raise FileNotFoundError(
            "Required tools not found in PATH: " + ", ".join(missing)
        )


def tool_version(tool_name: str) -> Optional[str]:
    """Best-effort first line of ``tool --version``, or None if unavailable."""
    try:
        result = subprocess.run(
            [tool_name, "--version"], capture_output=True, text=True, check=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    output = (result.stdout or "").strip() or (result.stderr or "").strip()
    if not output:
        return None
    return output.splitlines()[0].strip() or None


class ToolBase(ABC):
    def __init__(self, tool_name: str):
        self.tool_name = tool_name
        self.logger = get_logger(self.__class__.__name__)

    def _check_tool_installed(self) -> bool:
        if shutil.which(self.tool_name) is None:
            self.logger.error(f"{self.tool_name} not found in PATH")
            return False
        return True

    def _run_command(self, cmd: list, cwd: str = None) -> bool:
        try:
            self.logger.info(f"Running: {' '.join(cmd)}")
            subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=True)
            self.logger.info("Command completed successfully")
            return True
        except FileNotFoundError as e:
            raise RuntimeError(f"{self.tool_name} not found in PATH: {e}") from e
        except subprocess.CalledProcessError as e:
            stderr = (e.stderr or "").strip()
            self.logger.error(f"Command failed (exit {e.returncode}): {stderr}")
            raise RuntimeError(
                f"{self.tool_name} failed with exit code {e.returncode}: "
                f"{' '.join(map(str, e.cmd))}\n{stderr}"
            ) from e
    
    @abstractmethod
    def validate_inputs(self) -> bool:
        pass


class Tool(ToolBase):
    """Concrete implementation of ToolBase"""
    def validate_inputs(self) -> bool:
        return self._check_tool_installed()
