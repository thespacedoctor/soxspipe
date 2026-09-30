"""Deterministic, isolated defaults for the required test suite."""

from __future__ import annotations

import os
import socket
import tempfile
from pathlib import Path
from typing import NoReturn

import pytest

_SESSION_RUNTIME = tempfile.TemporaryDirectory(prefix="soxspipe-tests-")
_SESSION_ROOT = Path(_SESSION_RUNTIME.name)
_CATEGORY_MARKERS = frozenset({"unit", "integration", "real_data"})
_ORIGINAL_CREATE_CONNECTION = socket.create_connection
_ORIGINAL_SOCKET_CONNECT = socket.socket.connect
_ORIGINAL_SOCKET_CONNECT_EX = socket.socket.connect_ex
_ORIGINAL_SOCKET_SENDTO = socket.socket.sendto
_ORIGINAL_SOCKET_SENDMSG = getattr(socket.socket, "sendmsg", None)
_RESOLVER_NAMES = (
    "getaddrinfo",
    "getfqdn",
    "gethostbyaddr",
    "gethostbyname",
    "gethostbyname_ex",
    "getnameinfo",
)
_ORIGINAL_RESOLVERS = {
    resolverName: getattr(socket, resolverName)
    for resolverName in _RESOLVER_NAMES
    if hasattr(socket, resolverName)
}
_ENVIRONMENT_OVERRIDES = {
    "BLAS_NUM_THREADS": "1",
    "HOME": str(_SESSION_ROOT / "home"),
    "MKL_NUM_THREADS": "1",
    "MPLBACKEND": "Agg",
    "MPLCONFIGDIR": str(_SESSION_ROOT / "matplotlib"),
    "NUMBA_CACHE_DIR": str(_SESSION_ROOT / "numba"),
    "NUMEXPR_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "PYTHONDONTWRITEBYTECODE": "1",
    "TERM": "vt100",
    "VECLIB_MAXIMUM_THREADS": "1",
    "XDG_CACHE_HOME": str(_SESSION_ROOT / "cache"),
    "XDG_CONFIG_HOME": str(_SESSION_ROOT / "config"),
}
_ORIGINAL_ENVIRONMENT = {
    variableName: os.environ.get(variableName)
    for variableName in _ENVIRONMENT_OVERRIDES
}


def _deny_network(*args: object, **kwargs: object) -> NoReturn:
    raise RuntimeError("Network access is disabled in required tests")


def _guarded_connect(socketObject: socket.SocketType, address: object) -> None:
    if socketObject.family in {socket.AF_INET, socket.AF_INET6}:
        _deny_network()
    _ORIGINAL_SOCKET_CONNECT(socketObject, address)


def _guarded_connect_ex(socketObject: socket.SocketType, address: object) -> int:
    if socketObject.family in {socket.AF_INET, socket.AF_INET6}:
        _deny_network()
    return _ORIGINAL_SOCKET_CONNECT_EX(socketObject, address)


def _guarded_sendto(
    socketObject: socket.SocketType,
    data: bytes,
    *args: object,
) -> int:
    if socketObject.family in {socket.AF_INET, socket.AF_INET6}:
        _deny_network()
    return _ORIGINAL_SOCKET_SENDTO(socketObject, data, *args)


def _guarded_sendmsg(
    socketObject: socket.SocketType,
    buffers: object,
    ancdata: object = (),
    flags: int = 0,
    address: object | None = None,
) -> int:
    if socketObject.family in {socket.AF_INET, socket.AF_INET6}:
        _deny_network()
    if _ORIGINAL_SOCKET_SENDMSG is None:
        raise AttributeError("socket.sendmsg is unavailable")
    if address is None:
        return _ORIGINAL_SOCKET_SENDMSG(socketObject, buffers, ancdata, flags)
    return _ORIGINAL_SOCKET_SENDMSG(socketObject, buffers, ancdata, flags, address)


for variableName, value in _ENVIRONMENT_OVERRIDES.items():
    os.environ[variableName] = value

socket.create_connection = _deny_network
for resolverName in _ORIGINAL_RESOLVERS:
    setattr(socket, resolverName, _deny_network)
socket.socket.connect = _guarded_connect
socket.socket.connect_ex = _guarded_connect_ex
socket.socket.sendto = _guarded_sendto
if _ORIGINAL_SOCKET_SENDMSG is not None:
    socket.socket.sendmsg = _guarded_sendmsg


class RecordingLogger:
    """Small logger compatible with the interfaces exercised by unit tests."""

    def __init__(self) -> None:
        self.messages: list[tuple[str, str]] = []
        self.handlers: list[object] = []

    def _record(self, level: str, message: object) -> None:
        self.messages.append((level, str(message)))

    def debug(self, message: object) -> None:
        self._record("debug", message)

    def info(self, message: object) -> None:
        self._record("info", message)

    def warning(self, message: object) -> None:
        self._record("warning", message)

    def error(self, message: object) -> None:
        self._record("error", message)

    def critical(self, message: object) -> None:
        self._record("critical", message)

    def print(self, message: object) -> None:
        self._record("print", message)


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Require one execution category on every collected test."""
    for item in items:
        categories = {
            marker.name
            for marker in item.iter_markers()
            if marker.name in _CATEGORY_MARKERS
        }
        if len(categories) != 1:
            raise pytest.UsageError(
                f"{item.nodeid} must declare exactly one execution category"
            )


def pytest_unconfigure() -> None:
    """Restore process-wide state and release the session runtime directory."""
    socket.create_connection = _ORIGINAL_CREATE_CONNECTION
    for resolverName, originalResolver in _ORIGINAL_RESOLVERS.items():
        setattr(socket, resolverName, originalResolver)
    socket.socket.connect = _ORIGINAL_SOCKET_CONNECT
    socket.socket.connect_ex = _ORIGINAL_SOCKET_CONNECT_EX
    socket.socket.sendto = _ORIGINAL_SOCKET_SENDTO
    if _ORIGINAL_SOCKET_SENDMSG is not None:
        socket.socket.sendmsg = _ORIGINAL_SOCKET_SENDMSG
    for variableName, originalValue in _ORIGINAL_ENVIRONMENT.items():
        if originalValue is None:
            os.environ.pop(variableName, None)
        else:
            os.environ[variableName] = originalValue
    _SESSION_RUNTIME.cleanup()


@pytest.fixture(autouse=True)
def isolated_runtime(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep runtime state out of the real home and package directories."""
    homePath = tmp_path / "home"
    configPath = tmp_path / "config"
    cachePath = tmp_path / "cache"
    matplotlibPath = cachePath / "matplotlib"
    numbaPath = cachePath / "numba"
    for path in (homePath, configPath, cachePath, matplotlibPath, numbaPath):
        path.mkdir(parents=True, exist_ok=True)

    monkeypatch.setenv("HOME", str(homePath))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(configPath))
    monkeypatch.setenv("XDG_CACHE_HOME", str(cachePath))
    monkeypatch.setenv("MPLCONFIGDIR", str(matplotlibPath))
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(numbaPath))
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")
    monkeypatch.chdir(tmp_path)


@pytest.fixture(autouse=True)
def block_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Reject accidental network access from required tests."""
    monkeypatch.setattr(socket, "create_connection", _deny_network)
    for resolverName in _ORIGINAL_RESOLVERS:
        monkeypatch.setattr(socket, resolverName, _deny_network)
    monkeypatch.setattr(socket.socket, "connect", _guarded_connect)
    monkeypatch.setattr(socket.socket, "connect_ex", _guarded_connect_ex)
    monkeypatch.setattr(socket.socket, "sendto", _guarded_sendto)
    if _ORIGINAL_SOCKET_SENDMSG is not None:
        monkeypatch.setattr(socket.socket, "sendmsg", _guarded_sendmsg)


@pytest.fixture
def log() -> RecordingLogger:
    """Return a fresh in-memory logger."""
    return RecordingLogger()
