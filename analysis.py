#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Data Analysis Package of 3D Relativistic Electromagnetic Hydrodyanmics

 $Id: analysis.py,v c3b95b77a55e 2014/09/12 09:02:19 amano $
"""

import os
import sys
import numpy as np
import ioutils

def debug_read(fp):
    rank = fp.read(np.int32)
    dims = fp.read(np.int32, (rank,))
    data = fp.read(np.float64, dims)
    return data

def stream_function(bx, by, delx, dely):
    if not (bx.shape == by.shape):
        raise RuntimeError
    Ny, Nx = bx.shape
    Iy, Ix = np.mgrid[1:Ny,1:Nx]
    phi = np.cumsum(bx*dely, axis=0) - np.cumsum(by*delx, axis=1)
    return 0.25*(phi[Iy,Ix] + phi[Iy-1,Ix] + phi[Iy,Ix-1] + phi[Iy-1,Ix-1])

class DataReader:
    "Data reader object"
    def __init__(self, filename):
        self.fp = ioutils.BinaryReader(filename)
        # read number of record
        self.fp.seek(-4, 2)
        self.Nt = self.fp.read(np.int32)
        # start from the head
        self.fp.seek(0, 0)
        self.Nbc = self.fp.read(np.int32)
        self.Nx  = self.fp.read(np.int32)
        self.Ny  = self.fp.read(np.int32)
        self.Nz  = self.fp.read(np.int32)
        self.delx = self.fp.read(np.float64)
        self.dely = self.fp.read(np.float64)
        self.delz = self.fp.read(np.float64)
        self.xig  = self.fp.read(np.float64, (self.Nx,))
        self.yig  = self.fp.read(np.float64, (self.Ny,))
        self.zig  = self.fp.read(np.float64, (self.Nz,))
        # physical parameters
        self.cc   = self.fp.read(np.float64)
        self.gam  = self.fp.read(np.float64)
        # proton
        self.qp   = self.fp.read(np.float64)
        self.mp   = self.fp.read(np.float64)
        # electron
        self.qe   = self.fp.read(np.float64)
        self.me   = self.fp.read(np.float64)
        #
        self.time = np.zeros((self.Nt,), dtype=np.float64)
        self.delt = np.zeros((self.Nt,), dtype=np.float64)
        dims1 = (self.Nz,self.Ny,self.Nx,10)
        dims2 = (self.Nz,self.Ny,self.Nx,6)
        dims3 = (self.Nz+1,self.Ny+1,self.Nx+1,6)
        self.uf   = np.zeros((self.Nt,) + dims1, dtype=np.float64)
        self.ebc  = np.zeros((self.Nt,) + dims2, dtype=np.float64)
        self.ueb  = np.zeros((self.Nt,) + dims3, dtype=np.float64)
        for n in range(self.Nt):
            self.time[n] = self.fp.read(np.float64)
            self.delt[n] = self.fp.read(np.float64)
            self.uf[n]   = self.fp.read(np.float64, dims1)
            self.ebc[n]  = self.fp.read(np.float64, dims2)
            self.ueb[n]  = self.fp.read(np.float64, dims3)
            count        = self.fp.read(np.int32)
        # alias
        self.np  = self.uf[...,0]
        self.upx = self.uf[...,1]
        self.upy = self.uf[...,2]
        self.upz = self.uf[...,3]
        self.pp  = self.uf[...,4]
        self.ne  = self.uf[...,5]
        self.uex = self.uf[...,6]
        self.uey = self.uf[...,7]
        self.uez = self.uf[...,8]
        self.pe  = self.uf[...,9]
        self.ex  = self.ueb[...,0]
        self.ey  = self.ueb[...,1]
        self.ez  = self.ueb[...,2]
        self.bx  = self.ueb[...,3]
        self.by  = self.ueb[...,4]
        self.bz  = self.ueb[...,5]

    def divE(self, It):
        Iz, Iy, Ix = np.mgrid[1:self.Nz+1,1:self.Ny+1,1:self.Nx+1]
        Ex = self.ex[It]
        Ey = self.ey[It]
        Ez = self.ez[It]
        gp = np.sqrt(1 + self.upx[It]**2 + self.upy[It]**2 + self.upz[It]**2)
        ge = np.sqrt(1 + self.uex[It]**2 + self.uey[It]**2 + self.uez[It]**2)
        ro = 4*np.pi*(self.qp*gp*self.np[It] + self.qe*ge*self.ne[It])
        err =\
            (Ex[Iz,Iy,Ix] - Ex[Iz,Iy,Ix-1])/self.delx + \
            (Ey[Iz,Iy,Ix] - Ey[Iz,Iy-1,Ix])/self.dely + \
            (Ez[Iz,Iy,Ix] - Ez[Iz-1,Iy,Ix])/self.delz - ro
        return err

    def divB(self, It):
        Iz, Iy, Ix = np.mgrid[1:self.Nz+1,1:self.Ny+1,1:self.Nx+1]
        Bx = self.bx[It]
        By = self.by[It]
        Bz = self.bz[It]
        err =\
            (Bx[Iz,Iy,Ix] - Bx[Iz,Iy,Ix-1])/self.delx + \
            (By[Iz,Iy,Ix] - By[Iz,Iy-1,Ix])/self.dely + \
            (Bz[Iz,Iy,Ix] - Bz[Iz-1,Iy,Ix])/self.delz
        return err

class DebugReader:
    "Debug data reader object"
    def __init__(self, filename):
        self.fp = ioutils.BinaryReader(filename)
        # read number of record
        self.fp.seek(-4, 2)
        self.Nt = self.fp.read(np.int32)
        # start from the head
        self.fp.seek(0, 0)
        self.Nbc = self.fp.read(np.int32)
        self.Nx  = self.fp.read(np.int32)
        self.Ny  = self.fp.read(np.int32)
        self.Nz  = self.fp.read(np.int32)
        self.delx = self.fp.read(np.float64)
        self.dely = self.fp.read(np.float64)
        self.delz = self.fp.read(np.float64)
        self.xig  = self.fp.read(np.float64, (self.Nx,))
        self.yig  = self.fp.read(np.float64, (self.Ny,))
        self.zig  = self.fp.read(np.float64, (self.Nz,))
        # physical parameters
        self.cc   = self.fp.read(np.float64)
        self.gam  = self.fp.read(np.float64)
        # proton
        self.qp   = self.fp.read(np.float64)
        self.mp   = self.fp.read(np.float64)
        # electron
        self.qe   = self.fp.read(np.float64)
        self.me   = self.fp.read(np.float64)

    def read_next(self):
        self.time  = self.fp.read(np.float64)
        self.delt  = self.fp.read(np.float64)
        self.uf    = debug_read(self.fp)
        self.ff    = debug_read(self.fp)
        self.ebc   = debug_read(self.fp)
        self.ueb   = debug_read(self.fp)
        self.feb   = debug_read(self.fp)
        self.div   = debug_read(self.fp)
        count = self.fp.read(np.int32)
