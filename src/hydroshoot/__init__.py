from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("hydroshoot")
except PackageNotFoundError:
    # package is not installed
    pass