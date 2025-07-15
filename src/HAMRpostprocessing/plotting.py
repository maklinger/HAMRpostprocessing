# functions to help with plotting
import numpy as np

# def plc_new(myvar, xcoord=None, ycoord=None, ax=None, **kwargs):  # plc
#     global r, h, ph
#     l = [None] * nb2d

#     if (np.min(myvar) == np.max(myvar)):
#         print("The quantity you are trying to plot is a constant = %g." % np.min(myvar))
#         return
#     cb = kwargs.pop('cb', False)
#     nc = kwargs.pop('nc', 15)
#     k = kwargs.pop('k', 0)
#     mirrory = kwargs.pop('mirrory', 0)
#     # cmap = kwargs.pop('cmap',cm.jet)
#     isfilled = kwargs.pop('isfilled', False)
#     xy = kwargs.pop('xy', 0)
#     xmax = kwargs.pop('xmax', 10)
#     ymax = kwargs.pop('ymax', 5)
#     z = kwargs.pop('z', 0)

#     if ax is None:
#         ax = plt.gca()
#     if isfilled:
#         for i in range(0, nb):
#             index_z_block=int((z-int((z/360))*360.0)/360.0*bs3new*nb3*(1+REF_3)**(block[n_ord[i], AMR_LEVEL3]))
#             if (block[n_ord[i], AMR_COORD3] == int(index_z_block/bs3new)):
#                 offset=index_z_block-block[n_ord[i], AMR_COORD3]*bs3new
#                 res = ax.contourf(xcoord[i, :, :, offset], ycoord[i, :, :, offset], myvar[i, :, :, offset], nc, extend='both',**kwargs)
#     else:
#         for i in range(0, nb):
#             index_z_block=int(z/360.0*bs3new*nb3*(1+REF_3)**(block[n_ord[i], AMR_LEVEL3]))
#             if (block[n_ord[i], AMR_COORD3] == int(index_z_block/bs3new)):
#                 offset=index_z_block-block[n_ord[i], AMR_COORD3]*bs3new
#                 res = ax.contour(xcoord[i, :, :, offset], ycoord[i, :, :, offset], myvar[i, :, :, offset], nc, linewidths=4, extend='both', **kwargs)
#     if (cb == True):  # use color bar
#         plt.colorbar(res, ax=ax)
#     if xy:
#         plt.xlim(-xmax, xmax)
#         plt.ylim(-ymax, ymax)
#     return res



# def plc_cart(var, min, max, rmax, offset, name, label):
#     global aphi, r, h, ph, print_fieldlines,notebook, do_box, print_jetsheath, print_jetEHT, print_alfven
#     global beta_gamma_cut, beta_gamma_sq
#     global rho_cf, sigma_cf, Be_cf
#     fig = plt.figure(figsize=(64, 32))

#     X = r*np.sin(h)
#     Y = r*np.cos(h)
#     if(nb==1 and do_box==0):
#         X[:,:,0]=0.0*X[:,:,0]
#         X[:,:,bs2new-1]=0.0*X[:,:,bs2new-1]

#     plotmax = int(10*rmax * np.sqrt(2))

#     ilim = len(r[0, :, 0, 0]) - 1
#     for i in range(len(r[0, :, 0, 0])):
#         if r[0, i, 0, 0] > np.sqrt(2)*plotmax:
#             ilim = i
#             break

#     plt.subplot(1, 2, 1)
#     plc_new(np.log10((var))[:, 0:ilim], levels=np.arange(min, max, (max-min)/300.0), cb=0, isfilled=1, xcoord=X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=rmax, ymax=rmax)
#     res = plc_new(np.log10((var))[:, 0:ilim], levels=np.arange(min, max, (max-min)/300.0), cb=0, isfilled=1, xcoord=-1.0 * X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=rmax, ymax=rmax)
    
#     if (print_fieldlines == 1):
#         B0, BR, Bphi, Bz = convert123Rpz((B)[:, 0, :, :, :])
#         iBR = reinterp(BR, (-rmax, rmax, -rmax, rmax), bs1new, domirror=1, isasymmetric=1, method="cubic", domask=0)
#         iBz = reinterp(Bz, (-rmax, rmax, -rmax, rmax), bs1new, domirror=1, method="cubic", domask=0)
#         plane_mag = iBR**2+iBz**2
#         x = np.linspace(-rmax, rmax, bs1new)
#         z = np.linspace(-rmax, rmax, bs1new)
#         Xm, Ym = np.meshgrid(x, z)
#         plt.streamplot(Xm, Ym, iBR, iBz, color="k", arrowsize=2, density=[1,4], linewidth=2, broken_streamlines=False)        
#         # plc_new(aphi[:, 0:ilim], levels=np.arange(aphi[:, 0:ilim].min(), aphi[:, 0:ilim].max(), (aphi[:, 0:ilim].max()-aphi[:, 0:ilim].min())/20.0), cb=0,colors="black", isfilled=0, xcoord=X[:, 0:ilim], ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=rmax, ymax=rmax)
#         # plc_new(aphi[:, 0:ilim], levels=np.arange(aphi[:, 0:ilim].min(), aphi[:, 0:ilim].max(), (aphi[:, 0:ilim].max()-aphi[:, 0:ilim].min())/20.0), cb=0,colors="black", isfilled=0, xcoord=-1.0 * X[:, 0:ilim], ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=rmax, ymax=rmax)
    
#     if (print_jetsheath == 1):
#         sigma_plot_max = np.log10(sigma_cf + 0.0001)
#         sigma_plot_min = np.log10(sigma_cf - 0.0001)
#         plc_new(np.log10((bsq/rho))[:, 0:ilim], levels=np.arange(sigma_plot_min, sigma_plot_max, (sigma_plot_max-sigma_plot_min)/1.0), cb=0, colors="black", isfilled=0, xcoord=X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=rmax, ymax=rmax)
#         plc_new(np.log10((bsq/rho))[:, 0:ilim], levels=np.arange(sigma_plot_min, sigma_plot_max, (sigma_plot_max-sigma_plot_min)/1.0), cb=0, colors="black",  isfilled=0, xcoord=-1 * X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=rmax, ymax=rmax)
    
#         Be_plot_max = np.log10(Be_cf + 0.0001)
#         Be_plot_min = np.log10(Be_cf - 0.0001)
#         plc_new(np.log10(np.abs(Be))[:, 0:ilim], levels=np.arange(Be_plot_min, Be_plot_max, (Be_plot_max-Be_plot_min)/1.0), cb=0, colors="magenta", isfilled=0, xcoord=X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=rmax, ymax=rmax)
#         plc_new(np.log10(np.abs(Be))[:, 0:ilim], levels=np.arange(Be_plot_min, Be_plot_max, (Be_plot_max-Be_plot_min)/1.0), cb=0, colors="magenta", isfilled=0, xcoord=-1 * X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=rmax, ymax=rmax)          
        
#     if (print_alfven == 1):
#         alfven_plot_max = np.log10(1.02 + 1e-4)
#         alfven_plot_min = np.log10(1.02 - 1e-4)
#         plc_new(np.log10(np.linalg.norm(uu[1:,0,:,:,:]*np.sqrt(gcov[1,1])/uu[0,0,:,:,:],axis=0)/((B[2]*np.sqrt(gcov[2,2])/(rho+(gam-1)*ug+bsq))**0.5))[:, 0:ilim], levels=np.arange(alfven_plot_min, alfven_plot_max, (alfven_plot_max-alfven_plot_min)/1.0), cb=0, colors="purple", isfilled=0, xcoord=X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=rmax, ymax=rmax)
#         plc_new(np.log10(np.linalg.norm(uu[1:,0,:,:,:]*np.sqrt(gcov[1,1])/uu[0,0,:,:,:],axis=0)/((B[2]*np.sqrt(gcov[2,2])/(rho+(gam-1)*ug+bsq))**0.5))[:, 0:ilim], levels=np.arange(alfven_plot_min, alfven_plot_max, (alfven_plot_max-alfven_plot_min)/1.0), cb=0, colors="purple",  isfilled=0, xcoord=-1 * X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=rmax, ymax=rmax)


#     plt.xlabel(r"$x / R_g$", fontsize=90)
#     plt.ylabel(r"$z / R_g$", fontsize=90)
#     plt.title(label, fontsize=90)
#     ax = plt.gca()
#     ax.yaxis.set_ticks_position('left')
#     ax.xaxis.set_ticks_position('bottom')
#     ax.tick_params(axis='both', reset=False, which='both', length=24, width=6)
#     plt.gca().set_aspect(1)
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.05)
#     cb=plt.colorbar(res, cax=cax)
#     #cb.ax.tick_params(labelsize=50)

#     plt.subplot(1, 2, 2)
#     plc_new(np.log10((var))[:, 0:ilim], levels=np.arange(min, max, (max-min)/300.0), cb=0, isfilled=1, xcoord=X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=rmax * 5, ymax=rmax * 5)
#     res = plc_new(np.log10((var))[:, 0:ilim], levels=np.arange(min, max, (max-min)/300.0), cb=0, isfilled=1, xcoord=-1.0 * X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=rmax * 5, ymax=rmax * 5)
#     if (print_fieldlines == 1):
#         B0, BR, Bphi, Bz = convert123Rpz((B)[:, 0, :, :, :])
#         iBR = reinterp(BR, (-5*rmax, 5*rmax, -5*rmax, 5*rmax), bs1new, domirror=1, isasymmetric=1, method="cubic", domask=0)
#         iBz = reinterp(Bz, (-5*rmax, 5*rmax, -5*rmax, 5*rmax), bs1new, domirror=1, method="cubic", domask=0)
#         x = np.linspace(-5*rmax, 5*rmax, bs1new)
#         z = np.linspace(-5*rmax, 5*rmax, bs1new)
#         Xm, Ym = np.meshgrid(x, z)
#         plt.streamplot(Xm, Ym, iBR, iBz, color="k", arrowsize=2, density=[1,4], linewidth=2, broken_streamlines=False)
        
#         # plc_new(aphi[:, 0:ilim], levels=np.arange(aphi[:, 0:ilim].min(), aphi[:, 0:ilim].max(), (aphi[:, 0:ilim].max()-aphi[:, 0:ilim].min())/20.0), cb=0,colors="black", isfilled=0, xcoord=X[:, 0:ilim], ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=rmax * 5, ymax=rmax * 5)
#         # plc_new(aphi[:, 0:ilim], levels=np.arange(aphi[:, 0:ilim].min(), aphi[:, 0:ilim].max(), (aphi[:, 0:ilim].max()-aphi[:, 0:ilim].min())/20.0), cb=0,colors="black", isfilled=0, xcoord=-1.0 * X[:, 0:ilim], ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=rmax * 5, ymax=rmax * 5)
#     if (print_jetsheath == 1):
#         sigma_plot_max = np.log10(sigma_cf + 0.0001)
#         sigma_plot_min = np.log10(sigma_cf - 0.0001)
#         plc_new(np.log10((bsq/rho))[:, 0:ilim], levels=np.arange(sigma_plot_min, sigma_plot_max, (sigma_plot_max-sigma_plot_min)/1.0), cb=0, colors="black", isfilled=0, xcoord=X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=rmax*5, ymax=rmax*5)
#         plc_new(np.log10((bsq/rho))[:, 0:ilim], levels=np.arange(sigma_plot_min, sigma_plot_max, (sigma_plot_max-sigma_plot_min)/1.0), cb=0, colors="black",  isfilled=0, xcoord=-1*X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=rmax*5, ymax=rmax*5)
    
#         Be_plot_max = np.log10(Be_cf + 0.0001)
#         Be_plot_min = np.log10(Be_cf - 0.0001)
#         plc_new(np.log10(np.abs(Be))[:, 0:ilim], levels=np.arange(Be_plot_min, Be_plot_max, (Be_plot_max-Be_plot_min)/1.0), cb=0, colors="magenta", isfilled=0, xcoord=X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=rmax*5, ymax=rmax*5)
#         plc_new(np.log10(np.abs(Be))[:, 0:ilim], levels=np.arange(Be_plot_min, Be_plot_max, (Be_plot_max-Be_plot_min)/1.0), cb=0, colors="magenta", isfilled=0, xcoord=-1*X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=rmax*5, ymax=rmax*5)
        
#     if (print_alfven == 1):
#         alfven_plot_max = np.log10(1. + 1e-4)
#         alfven_plot_min = np.log10(1. - 1e-4)
#         plc_new(np.log10(np.linalg.norm(uu[1:,0,:,:,:]*np.sqrt(gcov[1,1])/uu[0,0,:,:,:],axis=0)/((B[2]*np.sqrt(gcov[2,2])/(rho+(gam-1)*ug+bsq))**0.5))[:, 0:ilim], levels=np.arange(alfven_plot_min, alfven_plot_max, (alfven_plot_max-alfven_plot_min)/1.0), cb=0, colors="purple", isfilled=0, xcoord=X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=offset, xmax=5*rmax, ymax=5*rmax)
#         plc_new(np.log10(np.linalg.norm(uu[1:,0,:,:,:]*np.sqrt(gcov[1,1])/uu[0,0,:,:,:],axis=0)/((B[2]*np.sqrt(gcov[2,2])/(rho+(gam-1)*ug+bsq))**0.5))[:, 0:ilim], levels=np.arange(alfven_plot_min, alfven_plot_max, (alfven_plot_max-alfven_plot_min)/1.0), cb=0, colors="purple",  isfilled=0, xcoord=-1 * X[:, 0:ilim],ycoord=Y[:, 0:ilim], xy=1, z=180 + offset, xmax=5*rmax, ymax=5*rmax)

#     plt.xlabel(r"$x / R_g$", fontsize=90)
#     #plt.ylabel(r"$z / R_g$", fontsize=60)
#     plt.title(label, fontsize=90)
#     ax = plt.gca()
#     ax.yaxis.set_ticks_position('left')
#     ax.xaxis.set_ticks_position('bottom')
#     ax.tick_params(axis='both', reset=False, which='both', length=24, width=6)
#     plt.gca().set_aspect(1)
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.05)
#     cb=plt.colorbar(res, cax=cax)
#     #cb.ax.tick_params(labelsize=50)
#     plt.savefig(name, dpi=100)
#     if (notebook==0):
#         plt.close('all')
