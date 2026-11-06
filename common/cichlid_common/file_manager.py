"""The lab's FileManager, found where it is checked out. The real implementation is still helper_modules/file_manager.py of the original layout (the alignment run on the server uses an older copy of it, which can
overwrite another machine's sheet changes, so it is NOT moved until that run is finished: see MIGRATE.md). Every module imports it from here, so the move changes one file.
    from cichlid_common.file_manager import FileManager"""
import importlib
from .common_light import ensure_lab_modules_on_path


def _impl():
    ensure_lab_modules_on_path()
    last = None
    for name in ('cichlid_common.file_manager_impl', 'helper_modules.file_manager'):
        try:
            return importlib.import_module(name).FileManager
        except ImportError as e:
            last = e
    raise ImportError('FileManager not found (looked for cichlid_common.file_manager_impl and helper_modules.file_manager): ' + str(last))


def FileManager(*args, **kwargs):
    return _impl()(*args, **kwargs)


def FileManagerClass():
    return _impl()
