import shutil
import subprocess
from abc import ABC, abstractmethod
from rskit.utils.logger import get_logger

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
