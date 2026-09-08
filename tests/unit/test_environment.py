"""Tests for required-suite isolation."""

from __future__ import annotations

import inspect
import os
import socket
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit


def test_runtime_directories_are_isolated_under_tmp_path(tmp_path: Path) -> None:
    for variableName in (
        "HOME",
        "XDG_CONFIG_HOME",
        "XDG_CACHE_HOME",
        "MPLCONFIGDIR",
        "NUMBA_CACHE_DIR",
    ):
        assert Path(os.environ[variableName]).is_relative_to(tmp_path)


def test_numerical_libraries_are_limited_to_one_thread() -> None:
    for variableName in (
        "BLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    ):
        assert os.environ[variableName] == "1"


def test_network_connections_are_rejected() -> None:
    with pytest.raises(RuntimeError, match="Network access is disabled"):
        socket.create_connection(("example.com", 443))


def test_socket_class_remains_subclassable() -> None:
    assert inspect.isclass(socket.socket)

    class LateSocketSubclass(socket.socket):
        pass

    assert issubclass(LateSocketSubclass, socket.socket)


def test_name_resolution_is_rejected() -> None:
    with pytest.raises(RuntimeError, match="Network access is disabled"):
        socket.getaddrinfo("localhost", 443)


@pytest.mark.skipif(not hasattr(socket, "AF_UNIX"), reason="AF_UNIX unavailable")
def test_unix_socket_preserves_fileno_construction_semantics() -> None:
    originalSocket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    fileDescriptor = originalSocket.detach()
    wrappedSocket = None
    try:
        wrappedSocket = socket.socket(fileno=fileDescriptor)
        assert wrappedSocket.family == socket.AF_UNIX
        assert wrappedSocket.fileno() == fileDescriptor
    finally:
        if wrappedSocket is not None:
            wrappedSocket.close()
        else:
            os.close(fileDescriptor)


def test_udp_sendto_is_rejected() -> None:
    internetSocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        with pytest.raises(RuntimeError, match="Network access is disabled"):
            internetSocket.sendto(b"blocked", ("127.0.0.1", 9))
    finally:
        internetSocket.close()


@pytest.mark.skipif(
    not hasattr(socket.SocketType, "sendmsg"),
    reason="socket.sendmsg unavailable",
)
def test_udp_sendmsg_is_rejected() -> None:
    internetSocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        with pytest.raises(RuntimeError, match="Network access is disabled"):
            internetSocket.sendmsg([b"blocked"], [], 0, ("127.0.0.1", 9))
    finally:
        internetSocket.close()
