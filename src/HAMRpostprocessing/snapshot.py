import numpy as np
import os
import pathlib
import pp_c2

# try:
#     import pp_c
# except ImportError:
#     from c_ext.build_c_extensions import build_pp_c
#     build_pp_c()
#     import pp_c








class Snapshot:

    def __init__(self, path, dump_nr, 
                 axisym=1, interpolate_var=0, export_raytracing_RAZIEH=0, do_box=0,
                 lowres1=1, lowres2=1, lowres3=1, check_files=0, DISK_THICKNESS=0.02, verbosity=0):
        self.verbosity = verbosity
        self.path = pathlib.Path(path)
        self.dump_nr = dump_nr
        self.axisym = axisym
        self.interpolate_var = interpolate_var
        self.export_raytracing_RAZIEH = export_raytracing_RAZIEH
        self.do_box = do_box
        self.lowres1 = lowres1
        self.lowres2 = lowres2
        self.lowres3 = lowres3
        self.check_files = check_files
        self.DISK_THICKNESS = DISK_THICKNESS


        self.rank = 0
        self.mytype = np.float32
        self.AMR_ACTIVE = 0
        self.AMR_LEVEL = 1
        self.AMR_REFINED = 2
        self.AMR_COORD1 = 3
        self.AMR_COORD2 = 4
        self.AMR_COORD3 = 5
        self.AMR_PARENT = 6
        self.AMR_CHILD1 = 7
        self.AMR_CHILD2 = 8
        self.AMR_CHILD3 = 9
        self.AMR_CHILD4 = 10
        self.AMR_CHILD5 = 11
        self.AMR_CHILD6 = 12
        self.AMR_CHILD7 = 13
        self.AMR_CHILD8 = 14
        self.AMR_NBR1 = 15
        self.AMR_NBR2 = 16
        self.AMR_NBR3 = 17
        self.AMR_NBR4 = 18
        self.AMR_NBR5 = 19
        self.AMR_NBR6 = 20
        self.AMR_NODE = 21
        self.AMR_POLE = 22
        self.AMR_GROUP = 23
        self.AMR_CORN1 = 24
        self.AMR_CORN2 = 25
        self.AMR_CORN3 = 26
        self.AMR_CORN4 = 27
        self.AMR_CORN5 = 28
        self.AMR_CORN6 = 29
        self.AMR_CORN7 = 30
        self.AMR_CORN8 = 31
        self.AMR_CORN9 = 32
        self.AMR_CORN10 = 33
        self.AMR_CORN11 = 34
        self.AMR_CORN12 = 35
        self.AMR_LEVEL1=  110
        self.AMR_LEVEL2 = 111
        self.AMR_LEVEL3 = 112
        self.AMR_NBR1_3=113
        self.AMR_NBR1_4=114
        self.AMR_NBR1_7=115
        self.AMR_NBR1_8=116
        self.AMR_NBR2_1=117
        self.AMR_NBR2_2=118
        self.AMR_NBR2_3=119
        self.AMR_NBR2_4=120
        self.AMR_NBR3_1=121
        self.AMR_NBR3_2=122
        self.AMR_NBR3_5=123
        self.AMR_NBR3_6=124
        self.AMR_NBR4_5=125
        self.AMR_NBR4_6=126
        self.AMR_NBR4_7=127
        self.AMR_NBR4_8=128
        self.AMR_NBR5_1=129
        self.AMR_NBR5_3=130
        self.AMR_NBR5_5=131
        self.AMR_NBR5_7=132
        self.AMR_NBR6_2=133
        self.AMR_NBR6_4=134
        self.AMR_NBR6_6=135
        self.AMR_NBR6_8=136
        self.AMR_NBR1P=161
        self.AMR_NBR2P=162
        self.AMR_NBR3P=163
        self.AMR_NBR4P=164
        self.AMR_NBR5P=165
        self.AMR_NBR6P=166
        self.AMR_TIMELEVEL=36

    # def load_snapshot(self):
    #     self.read_block()
    #     self.read_parameters()
    #     self.read_grid_dumps()
    #     self.read_primitive_dumps()


    def read_block(self):
        # rblock_new
        # Read in data for every block
        if self.verbosity:
            print("read block...")
        oldgridfile = self.path / f"dumps{self.dump_nr:d}/grid"
        newgridfile = self.path / "gdumps/grid"
        if (os.path.isfile(oldgridfile)):
            self.gridfile = oldgridfile
            fin = open(oldgridfile, "rb")
            size = os.path.getsize(oldgridfile)
            self.nmax = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
            NV = 36
        elif(os.path.isfile(newgridfile)):
            self.gridfile = newgridfile
            fin = open(newgridfile, "rb")
            size = os.path.getsize(newgridfile)
            self.nmax = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
            NV = (size - 1) // self.nmax // 4
        else:
            print("Cannot find grid file in dump %d !" % self.dump_nr)

        # Allocate memory
        self.block = np.zeros((self.nmax, 200), dtype=np.int32, order='C')
        self.n_ord = np.zeros((self.nmax), dtype=np.int32, order='C')

        gd = np.fromfile(fin, dtype=np.int32, count=NV * self.nmax, sep='')
        gd = gd.reshape((NV, self.nmax), order='F').T
        self.block[:,0:NV] = gd
        if(NV<170):
            self.block[:, self.AMR_LEVEL1] = gd[:, self.AMR_LEVEL]
            self.block[:, self.AMR_LEVEL2] = gd[:, self.AMR_LEVEL]
            self.block[:, self.AMR_LEVEL3] = gd[:, self.AMR_LEVEL]

        i = 0
        if (os.path.isfile(oldgridfile)):
            for n in range(0, self.nmax):
                if self.block[n, self.AMR_ACTIVE] == 1:
                    self.n_ord[i] = n
                    i += 1

        fin.close()



    def read_parameters(self):
        # rpar_new
        if self.verbosity:
            print("read parameters...")
        parameterfile = self.path / f"dumps{self.dump_nr:d}/parameters"
        if (os.path.isfile(parameterfile)):
            fin = open(parameterfile, "rb")
            self.parameterfile = parameterfile
        else:
            print("Rpar error!")

        self.t = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.n_active = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.n_active_total = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.nstep = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.Dtd = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.Dtl = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.Dtr = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.dump_cnt = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.rdump_cnt = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.dt = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.failed = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]

        self.bs1 = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.bs2 = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.bs3 = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.nmax = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.nb1 = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.nb2 = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        self.nb3 = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]

        self.startx1 = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.startx2 = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.startx3 = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self._dx1 = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self._dx2 = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self._dx3 = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.tf = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.a = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.gam = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.cour = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.Rin = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.Rout = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.R0 = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        self.density_scale = np.fromfile(fin, dtype=np.float64, count=1, sep='')[0]
        for n in range(0,13):
            trash = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
        trash = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]

        if (trash >= 10000000):
            self.DOTRACER = 1
            trash = trash - 10000000
        else:
            self.DOTRACER = 0    
        if (trash >= 1000000):
            self.DO_NUCLEAR = 1
            trash = trash - 1000000
        else:
            self.DO_NUCLEAR = 0
        if (trash >= 100000):
            self.NEUTRINOS_M1 = 1
            trash = trash - 100000
        else:
            self.NEUTRINOS_M1 = 0
        if (trash >= 10000):
            self.DO_YE = 1
            trash = trash - 10000
        else:
            self.DO_YE = 0
        if(trash >= 1000):
            self.P_NUM=1
            trash=trash-1000
        else:
            self.P_NUM=0
        if (trash >= 100):
            self.TWO_T = 1
            # print(TWO_T)
            trash = trash - 100
        else:
            self.TWO_T = 0
        if (trash >= 10):
            self.RESISTIVE = 1
            trash = trash - 10
        else:
            self.RESISTIVE = 0
        if (trash >= 1):
            self.RAD_M1 = 1
            trash = trash - 1
        else:
            self.RAD_M1 = 0
        trash = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]

        #Set grid spacing
        self._dx1=(np.log(self.Rout)-np.log(self.Rin))/(self.bs1*self.nb1)
        self.fractheta=-self.startx2
        #fractheta = 1.0 - 2.0 / (bs2*nb2) * (bs3*nb3>2.0)
        self._dx2=2.0*self.fractheta/(self.bs2*self.nb2)
        self._dx3=2.0*np.pi/(self.bs3*self.nb3)

        self.nb = self.n_active_total
        self.rhor = 1 + (1 - self.a ** 2) ** 0.5

        self.NODE=np.copy(self.n_ord)
        self.TIMELEVEL=np.copy(self.n_ord)

        self.REF_1=1
        self.REF_2=1
        self.REF_3=1
        self.flag_restore = 0
        self.size = os.path.getsize(parameterfile)
        if(self.size>=66*4+3*self.n_active_total*4):
            n=0
            while n<self.n_active_total:
                self.n_ord[n]=np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
                self.TIMELEVEL[n] = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
                self.NODE[n] = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
                n=n+1
        elif(self.size >= 66 * 4 + 2 * self.n_active_total * 4):
            n = 0
            self.flag_restore=1
            while n < self.n_active_total:
                self.n_ord[n] = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
                self.TIMELEVEL[n] = np.fromfile(fin, dtype=np.int32, count=1, sep='')[0]
                n = n + 1

        if(
            self.export_raytracing_RAZIEH==1 and (
                self.bs1%self.lowres1!=0 or self.bs2%self.lowres2!=0 or 
                self.bs3%self.lowres3!=0 or (
                    (self.lowres1 & (self.lowres1-1) == 0) and 
                    self.lowres1 != 0)!=1 or (
                        (self.lowres2 & (self.lowres2-1) == 0) and 
                        self.lowres2 != 0)!=1 or (
                            (self.lowres3 & (self.lowres3-1) == 0) and 
                            self.lowres3 != 0)!=1)):
            if(self.rank==0):
                print("For raytracing block size needs to be divisable by lowres!")
        if(self.export_raytracing_RAZIEH==1 and self.interpolate_var==0):
            if (self.rank == 0):
                print("Warning: Variable interpolation is highly recommended for raytracing!")
        fin.close()

    def read_grid_dumps(self):
        # rgdump_griddata
        
        if self.verbosity:
            print("read grid data...")
        # self.set_cart=0
        # self.set_xc=0

        self.ACTIVE1 = np.max(self.block[self.n_ord, self.AMR_LEVEL1])*self.REF_1
        self.ACTIVE2 = np.max(self.block[self.n_ord, self.AMR_LEVEL2])*self.REF_2
        self.ACTIVE3 = np.max(self.block[self.n_ord, self.AMR_LEVEL3])*self.REF_3

        if ((int(self.nb1 * (1 + self.REF_1) ** self.ACTIVE1 * self.bs1) % self.lowres1) != 0 or 
            (int(self.nb2 * (1 + self.REF_2) ** self.ACTIVE2 * self.bs2) % self.lowres2) != 0 or 
            (int(self.nb3 * (1 + self.REF_3) ** self.ACTIVE3 * self.bs3) % self.lowres3) != 0):
            print("Incompatible lowres settings in rgdump_griddata")

        self.gridsizex1 = int(self.nb1 * (1 + self.REF_1) ** self.ACTIVE1 * self.bs1/self.lowres1)
        self.gridsizex2 = int(self.nb2 * (1 + self.REF_2) ** self.ACTIVE2 * self.bs2/self.lowres2)
        self.gridsizex3 = int(self.nb3 * (1 + self.REF_3) ** self.ACTIVE3 * self.bs3/self.lowres3)

        self._dx1 = self._dx1 * self.lowres1 * (1.0 / (1.0 + self.REF_1) ** self.ACTIVE1)
        self._dx2 = self._dx2 * self.lowres2 * (1.0 / (1.0 + self.REF_2) ** self.ACTIVE2)
        self._dx3 = self._dx3 * self.lowres3 * (1.0 / (1.0 + self.REF_3) ** self.ACTIVE3)

        #Calculate inner and outer boundaries of selection box after upscaling and downscaling; Assumes uniform grid x1=log(r) etc
        # if(do_box==1):
        #     i_min = max(np.int32((np.log(r_min)-(startx1+0.5*_dx1)) / _dx1) + 1, 0)
        #     i_max = min(np.int32((np.log(r_max)-(startx1+0.5*_dx1)) / _dx1) + 1, gridsizex1)
        #     j_min=max(np.int32(((2.0/np.pi*(theta_min)-1.0)-(startx2+0.5*_dx2))/_dx2) + 1,0)
        #     j_max=min(np.int32(((2.0/np.pi*(theta_max)-1.0)-(startx2+0.5*_dx2))/_dx2) + 1,gridsizex2)
        #     z_min=max(np.int32((phi_min-(startx3+0.5*_dx3))/_dx3) + 1,0)
        #     z_max=min(np.int32((phi_max-(startx3+0.5*_dx3))/_dx3) + 1,gridsizex3)

        #     gridsizex1 = i_max-i_min
        #     gridsizex2 = j_max-j_min
        #     gridsizex3 = z_max-z_min

        #     if((j_max<j_min or i_max<i_min or z_max<z_min) and rank==0):
        #         print("Bad box selection")
        # else:
        self.i_min=0
        self.i_max=self.gridsizex1
        self.j_min=0
        self.j_max=self.gridsizex2
        self.z_min=0
        self.z_max=self.gridsizex3

        self.nx = self.gridsizex1
        self.ny = self.gridsizex2
        self.nz = self.gridsizex3

        # Allocate memory
        self.x1 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        self.x2 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        self.x3 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        self.r = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        self.h = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        self.ph = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')

        if self.axisym:
            self.gcov = np.zeros((4, 4, 1, self.gridsizex1, self.gridsizex2, 1), dtype=self.mytype, order='C')
            self.gcon = np.zeros((4, 4, 1, self.gridsizex1, self.gridsizex2, 1), dtype=self.mytype, order='C')
            self.gdet = np.zeros((1, self.gridsizex1, self.gridsizex2, 1), dtype=self.mytype, order='C')
            self.dxdxp = np.zeros((4, 4, 1, self.gridsizex1, self.gridsizex2, 1), dtype=self.mytype, order='C')
        else:
            self.gcov = np.zeros((4, 4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.gcon = np.zeros((4, 4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.gdet = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.dxdxp = np.zeros((4, 4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')

        
        if(self.rank==0 and self.check_files==1):
            for n in range(0,self.n_active_total):
                gdumpfile_n = self.path / f"gdumps/gdump{self.n_ord[n]:d}"
                if(os.path.isfile(gdumpfile_n)==0 
                   or os.path.getsize(gdumpfile_n)!=(
                       9*self.bs1*self.bs2*self.bs3+(self.bs1*self.bs2*49)*(self.axisym)+
                       (self.bs1*self.bs2*self.bs3*49)*(self.axisym==0))*8):
                    print("Gdump file %d doesn't exist" %self.n_ord[n])
        print(self.path / f"gdumps/gdump{self.n_ord[0]:d}")
        self.size = os.path.getsize(self.path / f"gdumps/gdump{self.n_ord[0]:d}")
        if(self.size==58*self.bs3*self.bs2*self.bs1*8):
            flag=1
        else:
            flag=0
        print(str(self.path))
        pp_c2.rgdump_griddata(
            flag, self.interpolate_var, str(self.path), self.axisym, 
            self.n_ord, self.lowres1, self.lowres2, self.lowres3,
            self.nb, self.bs1, self.bs2, self.bs3, self.x1, self.x2, self.x3,
            self.r, self.h, self.ph, self.gcov, self.gcon, self.dxdxp,
            self.gdet, self.block, self.nb1, self.nb2, self.nb3, 
            self.REF_1, self.REF_2, self.REF_3, 
            np.max(self.block[self.n_ord, self.AMR_LEVEL1]), 
            np.max(self.block[self.n_ord, self.AMR_LEVEL2]), 
            np.max(self.block[self.n_ord, self.AMR_LEVEL3]), 
            self.startx1, self.startx2, self.startx3, 
            self._dx1, self._dx2, self._dx3, 
            self.export_raytracing_RAZIEH, self.i_min, self.i_max,
            self.j_min, self.j_max, self.z_min, self.z_max)


    def read_primitive_dumps(self):
        # rdump_griddata
        
        if self.verbosity:
            print("read primitive variables...")
        
        # Allocate memory
        self.rho = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        self.ug = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        self.uu = np.zeros((4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        self.B = np.zeros((4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        if(self.export_raytracing_RAZIEH):
            self.Rdot = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        else:
            self.Rdot = np.zeros((1, 1, 1, 1), dtype=self.mytype, order='C')
        self.bsq = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')

        if (self.DO_YE):
            self.Ye = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        else:
            self.Ye = self.ug

        if (self.DO_NUCLEAR):
            self.X_alpha = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.X_atm = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        else:
            self.X_alpha = self.ug
            self.X_atm = self.ug

        if (self.NEUTRINOS_M1):
            self.E_nu1 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.E_nu2 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.E_nu3 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.uu_nu1 = np.zeros((4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.uu_nu2 = np.zeros((4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.uu_nu3 = np.zeros((4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.neutrino_number1 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.neutrino_number2 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.neutrino_number3 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        else:
            self.E_nu1 = self.rho
            self.E_nu2 = self.rho
            self.E_nu3 = self.rho
            self.uu_nu1 = self.uu
            self.uu_nu2 = self.uu
            self.uu_nu3 = self.uu
            self.neutrino_number1 = self.rho
            self.neutrino_number2 = self.rho
            self.neutrino_number3 = self.rho

        if(self.RAD_M1):
            self.E_rad = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.uu_rad = np.zeros((4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        else:
            self.E_rad=np.copy(self.ug)
            self.uu_rad=np.copy(self.uu)

        if (self.DOTRACER):
            self.index_tr1 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.index_tr2 = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        else:
            self.index_tr1=self.rho
            self.index_tr2=self.rho

        if (self.RESISTIVE):
            self.E = np.zeros((4, 1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        else:
            self.E = self.B

        if (self.TWO_T):
            self.TE = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
            self.TI = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        else:
            self.TE = self.rho
            self.TI = self.rho

        if (self.P_NUM):
            self.photon_number = np.zeros((1, self.gridsizex1, self.gridsizex2, self.gridsizex3), dtype=self.mytype, order='C')
        else:
            self.photon_number = self.rho

        olddatadumpfile = self.path / f"dumps{self.dump_nr:d}/new_dump"
        if (os.path.isfile(olddatadumpfile)):
            flag = 1
        else:
            print("flag=0")
            if (self.rank == 0 and self.check_files == 10):
                for count in range(0, 5400):
                    if (os.path.isfile(self.path + "/dumps%d/new_dump%d" %(self.dump_nr, count))==0):
                        print("Dump file %d in folder %d doesn't exist" %(count,self.dump_nr))
            flag = 0


        pp_c2.rdump_griddata(
            flag, self.interpolate_var, 
            np.int32(self.RAD_M1),np.int32(self.RESISTIVE), np.int32(self.DOTRACER), 
            np.int32(self.TWO_T), np.int32(self.P_NUM), np.int32(self.DO_YE), 
            np.int32(self.DO_NUCLEAR), np.int32(self.NEUTRINOS_M1), str(self.path),
            self.dump_nr, self.n_active_total, self.lowres1, self.lowres2, self.lowres3, 
            self.nb, self.bs1, self.bs2, self.bs3, 
            self.rho, self.index_tr1, self.index_tr2, self.ug, self.uu, self.B, self.E, 
            self.E_rad, self.uu_rad, self.TE, self.TI, self.photon_number, self.Ye, 
            self.X_alpha, self.X_atm, self.E_nu1, self.E_nu2, self.E_nu3, 
            self.uu_nu1, self.uu_nu2, self.uu_nu3, self.neutrino_number1,
            self.neutrino_number2, self.neutrino_number3, self.gcov, self.gcon,
            self.axisym, self.n_ord, self.block, self.nb1, self.nb2, self.nb3,
            self.REF_1, self.REF_2, self.REF_3, 
            np.max(self.block[self.n_ord, self.AMR_LEVEL1]),
            np.max(self.block[self.n_ord, self.AMR_LEVEL2]), 
            np.max(self.block[self.n_ord, self.AMR_LEVEL3]),
            self.export_raytracing_RAZIEH, self.DISK_THICKNESS, 
            self.a, self.gam, self.Rdot, self.bsq, self.r, self.startx1,
            self.startx2, self.startx3, self._dx1, self._dx2, self._dx3,
            self.x1, self.x2, self.x3, self.i_min, self.i_max, self.j_min, self.j_max, 
            self.z_min, self.z_max)

        # pp_c.rdump_griddata(
        #     flag, self.interpolate_var, np.int32(self.RAD_M1), 
        #     np.int32(self.RESISTIVE), np.int32(self.TWO_T), np.int32(self.P_NUM),
        #     np.int32(self.DO_YE), np.int32(self.DO_NUCLEAR), np.int32(self.NEUTRINOS_M1), 
        #     str(self.path), self.dump_nr, self.n_active_total, 
        #     self.lowres1, self.lowres2, self.lowres3, self.nb,
        #     self.bs1, self.bs2, self.bs3, self.rho, self.ug, self.uu, self.B, self.E, 
        #     self.E_rad, self.uu_rad, self.TE, self.TI, self.photon_number, self.Ye, 
        #     self.X_alpha, self.X_atm, self.E_nu1, self.E_nu2, self.E_nu3, self.uu_nu1, 
        #     self.uu_nu2, self.uu_nu3, self.neutrino_number1, self.neutrino_number2, 
        #     self.neutrino_number3, self.gcov, self.gcon, self.axisym, self.n_ord,
        #     self.block, self.nb1, self.nb2, self.nb3, self.REF_1, self.REF_2, self.REF_3,
        #     np.max(self.block[self.n_ord, self.AMR_LEVEL1]),
        #     np.max(self.block[self.n_ord, self.AMR_LEVEL2]), 
        #     np.max(self.block[self.n_ord, self.AMR_LEVEL3]),
        #     self.export_raytracing_RAZIEH, self.DISK_THICKNESS, self.a, 
        #     self.gam, self.Rdot, self.bsq, self.r, 
        #     self.startx1, self.startx2, self.startx3, self._dx1, self._dx2, self._dx3,
        #     self.x1, self.x2, self.x3, self.i_min, self.i_max, 
        #     self.j_min, self.j_max, self.z_min, self.z_max)

        self.bs1new = self.gridsizex1
        self.bs2new = self.gridsizex2
        self.bs3new = self.gridsizex3

        if (self.do_box == 1):
            self.startx1 = self.startx1 + (self.i_min) * self._dx1
            self.startx2 = self.startx2 + (self.j_min) * self._dx2
            self.startx3 = self.startx3 + (self.z_min) * self._dx3

        self.nb2d = self.nb
        self.nb = 1
        self.nb1 = 1
        self.nb2 = 1
        self.nb3 = 1
        



