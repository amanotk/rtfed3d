#!/usr/bin/env python3

"""Data analysis helpers for rtfed3d custom binary outputs."""

from pathlib import Path

import numpy as np


class BinaryReader:
    """Small binary reader compatible with the old ioutils.BinaryReader API."""

    def __init__(self, filename):
        self.path = Path(filename)
        self._fp = self.path.open("rb")

    def close(self):
        self._fp.close()

    def seek(self, offset, whence=0):
        self._fp.seek(offset, whence)

    def read(self, dtype, shape=None):
        dtype = np.dtype(dtype)

        if shape is None:
            count = 1
        elif isinstance(shape, tuple):
            count = int(np.prod(shape, dtype=np.int64))
        else:
            count = int(shape)

        data = np.fromfile(self._fp, dtype=dtype, count=count)
        if data.size != count:
            raise EOFError(
                f"Unexpected end of file while reading {dtype} from {self.path}"
            )

        if shape is None:
            return data[0]

        if isinstance(shape, tuple):
            return data.reshape(shape)

        return data

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()


def debug_read(fp):
    rank = int(fp.read(np.int32))
    dims = tuple(int(v) for v in fp.read(np.int32, (rank,)))
    data = fp.read(np.float64, dims)
    return data


def stream_function(bx, by, delx, dely):
    if bx.shape != by.shape:
        raise RuntimeError("bx and by must have the same shape")
    ny, nx = bx.shape
    iy, ix = np.mgrid[1:ny, 1:nx]
    phi = np.cumsum(bx * dely, axis=0) - np.cumsum(by * delx, axis=1)
    return 0.25 * (
        phi[iy, ix] + phi[iy - 1, ix] + phi[iy, ix - 1] + phi[iy - 1, ix - 1]
    )


class DataReader:
    """Data reader object."""

    def __init__(self, filename):
        self.fp = BinaryReader(filename)
        self.fp.seek(-4, 2)
        self.Nt = int(self.fp.read(np.int32))
        self.fp.seek(0, 0)
        self.Nbc = int(self.fp.read(np.int32))
        self.Nx = int(self.fp.read(np.int32))
        self.Ny = int(self.fp.read(np.int32))
        self.Nz = int(self.fp.read(np.int32))
        self.delx = float(self.fp.read(np.float64))
        self.dely = float(self.fp.read(np.float64))
        self.delz = float(self.fp.read(np.float64))
        self.xig = self.fp.read(np.float64, (self.Nx,))
        self.yig = self.fp.read(np.float64, (self.Ny,))
        self.zig = self.fp.read(np.float64, (self.Nz,))
        self.cc = float(self.fp.read(np.float64))
        self.gam = float(self.fp.read(np.float64))
        self.qp = float(self.fp.read(np.float64))
        self.mp = float(self.fp.read(np.float64))
        self.qe = float(self.fp.read(np.float64))
        self.me = float(self.fp.read(np.float64))

        self.time = np.zeros((self.Nt,), dtype=np.float64)
        self.delt = np.zeros((self.Nt,), dtype=np.float64)
        dims1 = (self.Nz, self.Ny, self.Nx, 10)
        dims2 = (self.Nz, self.Ny, self.Nx, 6)
        dims3 = (self.Nz + 1, self.Ny + 1, self.Nx + 1, 6)
        self.uf = np.zeros((self.Nt,) + dims1, dtype=np.float64)
        self.ebc = np.zeros((self.Nt,) + dims2, dtype=np.float64)
        self.ueb = np.zeros((self.Nt,) + dims3, dtype=np.float64)

        for n in range(self.Nt):
            self.time[n] = self.fp.read(np.float64)
            self.delt[n] = self.fp.read(np.float64)
            self.uf[n] = self.fp.read(np.float64, dims1)
            self.ebc[n] = self.fp.read(np.float64, dims2)
            self.ueb[n] = self.fp.read(np.float64, dims3)
            self.fp.read(np.int32)

        self.np = self.uf[..., 0]
        self.upx = self.uf[..., 1]
        self.upy = self.uf[..., 2]
        self.upz = self.uf[..., 3]
        self.pp = self.uf[..., 4]
        self.ne = self.uf[..., 5]
        self.uex = self.uf[..., 6]
        self.uey = self.uf[..., 7]
        self.uez = self.uf[..., 8]
        self.pe = self.uf[..., 9]
        self.ex = self.ueb[..., 0]
        self.ey = self.ueb[..., 1]
        self.ez = self.ueb[..., 2]
        self.bx = self.ueb[..., 3]
        self.by = self.ueb[..., 4]
        self.bz = self.ueb[..., 5]

    def close(self):
        self.fp.close()

    def divE(self, it):
        iz, iy, ix = np.mgrid[1 : self.Nz + 1, 1 : self.Ny + 1, 1 : self.Nx + 1]
        ex = self.ex[it]
        ey = self.ey[it]
        ez = self.ez[it]
        gp = np.sqrt(1 + self.upx[it] ** 2 + self.upy[it] ** 2 + self.upz[it] ** 2)
        ge = np.sqrt(1 + self.uex[it] ** 2 + self.uey[it] ** 2 + self.uez[it] ** 2)
        ro = 4 * np.pi * (self.qp * gp * self.np[it] + self.qe * ge * self.ne[it])
        err = (
            (ex[iz, iy, ix] - ex[iz, iy, ix - 1]) / self.delx
            + (ey[iz, iy, ix] - ey[iz, iy - 1, ix]) / self.dely
            + (ez[iz, iy, ix] - ez[iz - 1, iy, ix]) / self.delz
            - ro
        )
        return err

    def divB(self, it):
        iz, iy, ix = np.mgrid[1 : self.Nz + 1, 1 : self.Ny + 1, 1 : self.Nx + 1]
        bx = self.bx[it]
        by = self.by[it]
        bz = self.bz[it]
        err = (
            (bx[iz, iy, ix] - bx[iz, iy, ix - 1]) / self.delx
            + (by[iz, iy, ix] - by[iz, iy - 1, ix]) / self.dely
            + (bz[iz, iy, ix] - bz[iz - 1, iy, ix]) / self.delz
        )
        return err


class DebugReader:
    """Debug data reader object."""

    def __init__(self, filename):
        self.fp = BinaryReader(filename)
        self.fp.seek(-4, 2)
        self.Nt = int(self.fp.read(np.int32))
        self.fp.seek(0, 0)
        self.Nbc = int(self.fp.read(np.int32))
        self.Nx = int(self.fp.read(np.int32))
        self.Ny = int(self.fp.read(np.int32))
        self.Nz = int(self.fp.read(np.int32))
        self.delx = float(self.fp.read(np.float64))
        self.dely = float(self.fp.read(np.float64))
        self.delz = float(self.fp.read(np.float64))
        self.xig = self.fp.read(np.float64, (self.Nx,))
        self.yig = self.fp.read(np.float64, (self.Ny,))
        self.zig = self.fp.read(np.float64, (self.Nz,))
        self.cc = float(self.fp.read(np.float64))
        self.gam = float(self.fp.read(np.float64))
        self.qp = float(self.fp.read(np.float64))
        self.mp = float(self.fp.read(np.float64))
        self.qe = float(self.fp.read(np.float64))
        self.me = float(self.fp.read(np.float64))

    def close(self):
        self.fp.close()

    def read_next(self):
        self.time = self.fp.read(np.float64)
        self.delt = self.fp.read(np.float64)
        self.uf = debug_read(self.fp)
        self.ff = debug_read(self.fp)
        self.ebc = debug_read(self.fp)
        self.ueb = debug_read(self.fp)
        self.feb = debug_read(self.fp)
        self.div = debug_read(self.fp)
        self.fp.read(np.int32)
