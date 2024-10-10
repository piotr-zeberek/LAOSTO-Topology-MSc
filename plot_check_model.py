# basic imports
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd

# ============ for latex fonts ============
from matplotlib import rc  # , font_manager

rc(
    "text.latex", preamble=r"\usepackage{lmodern} \usepackage{physics}"
)  # this helps use the plots in tex files
plt.rcParams.update({"font.size": 11})
plt.rcParams.update(
    {
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.8,
        "ytick.minor.width": 0.8,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "axes.linewidth": 0.8,
        "lines.linewidth": 0.8,
        "patch.linewidth": 0.8,
        "legend.fontsize": 8,
        "legend.title_fontsize": 9,
        "legend.fancybox": False,
        "legend.frameon": False,
        "legend.handlelength": 1.0,
        # "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.4,
        "figure.titlesize": 12,
        "figure.figsize": [3.54, 3.0],
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "mathtext.fontset": "cm",
        "text.usetex": True,
        "font.family": "Computer Modern Roman",
    }
)
# ==========================================

###################

data_dir = Path("build/data")
plot_dir = Path("plot")

if not Path.exists(plot_dir):
    Path.mkdir(plot_dir)

###################


def load_data(filename, data_dir=data_dir):
    data = np.loadtxt(data_dir / filename)
    return (col for col in data.T)


def prep_image_data(x, y, z):
    # x = np.unique(x)
    # y = np.unique(y)

    # X, Y = np.meshgrid(x, y)
    # Z = z.reshape(len(x), len(y))

    # return (np.transpose(Z), [min(x), max(x), min(y), max(y)])
    data = pd.DataFrame({"x": x, "y": y, "z": z})
    data_pivot = data.pivot_table(index="y", columns="x", values="z")
    extent = [min(x), max(x), min(y), max(y)]
    return data_pivot, extent


def set_labels(xlabel=r"$x$", ylabel=r"$y$", title=""):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def reset_plot_setting():
    plt.close("all")
    # plt.figure(figsize=figsize, dpi=dpi)


def save_plot(filename, bbox_inches="tight", plot_dir=plot_dir, dpi=300):
    plt.tight_layout()
    plt.savefig(plot_dir / filename, bbox_inches=bbox_inches, dpi=dpi)
    reset_plot_setting()


##########################


def get_vertical_slice(x, y, z, axis, offset=0.0):
    ux = np.unique(x)
    uy = np.unique(y)

    match axis:
        case "x":
            diff = np.abs(ux - offset)
            val = ux[np.argmin(diff)]

            return uy, z[x == val]

        case "y":
            diff = np.abs(uy - offset)
            val = uy[np.argmin(diff)]

            return ux, z[y == val]

        case _:
            raise ValueError("Invalid axis")


# # slice
# k, *energies = load_data("BS.dat")
# for band_idx, E in enumerate(energies):
#     plt.plot(k, E)

# # plt.xlim(-0.2, 0.2)
# plt.ylim(-5, 5)
# set_labels(r"$k_x\ (1/a)$ ", r"$E$ (meV)", r"$\mu = 0.0$ meV, $k_y = 0.0$")
# plt.grid()
# save_plot("BS.png")

# #2D
# kx, ky, *energies = load_data("BS.dat")
# ukx, uky = np.unique(kx), np.unique(ky)
# x,y = np.meshgrid(ukx, uky)

# fig = plt.figure(figsize=(7, 7))
# ax = fig.add_subplot(111, projection="3d")

# for i, e in enumerate(energies):
#     ax.plot_surface(x,y, e.reshape(len(ukx), len(uky)).T, alpha=0.5, label=f"{i}")

# set_labels(r"$k_x$", r"$k_y$", r"$E$ (meV)")
# ax.legend(title="Band index")
# # save_plot("bands_3d.png", reset=False)
# plt.show()

# fig, axes = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(5, 5.5))

# for i, mu in enumerate(np.linspace(-2.93, -2.85, 9)):
#     k, *energies = load_data(f"BS{i}.dat")

#     for band_idx, E in enumerate(energies):
#         axes.flat[i].plot(k, E)

#     axes.flat[i].set_title(f"$\mu = {mu:.2f}$ meV")

# for ax in axes.flat:
#     ax.set_ylim(-0.3, 0.3)
#     ax.set_xlabel(r"$k_x$ ($1/a$)")
#     ax.set_ylabel(r"$E$ (meV)")
#     ax.set_box_aspect(1)
#     ax.label_outer()
#     ax.grid()

# fig.suptitle(r"Transition from trivial to topological phase for $\mu \sim -2.89$ meV")
# save_plot(f"transition.png")

# Chern numbers vs mu
mu, *Cs = load_data("mu_Cs.dat")
# mu, *Cs = np.loadtxt("LAOSTO/mu_Cs/denser/mu_Cs.dat", unpack=True)

# for i in range(3):
#     C = np.sum(Cs[i * 2 : (i + 1) * 2], axis=0)
#     plt.plot(mu, C, label=i)

# for i,C in enumerate(Cs):
#     plt.plot(mu, C, label=f"{i}")

plt.plot(mu, np.sum(Cs[:2], axis=0))

set_labels(r"$\mu$ (meV)", r"$C$", r"Chern Number for BZ grid $500 \times 500$")
plt.grid(alpha=0.5)
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Band")
save_plot("mu_Cs.png")

# # Chern numbers vs mu link and BC
# fig, axes = plt.subplots(1, 2, sharex=True, figsize=(5, 2.7))

# for i, filename in enumerate(["mu_Cs.dat", "mu_Cs_BC.dat"]):
#     mu, *Cs = load_data(filename)
#     Cs = np.sum(Cs[:2], axis=0)
#     axes[i].plot(mu, Cs)
#     axes[i].set_xlabel(r"$\mu$ (meV)")
#     axes[i].set_ylabel(r"$C$")
#     axes[i].set_title("From link variable" if not i else "From Berry curvature")
#     axes[i].grid()
#     axes[i].set_box_aspect(1)

# fig.suptitle(r"Chern Number for BZ grid $1500 \times 1500$")
# save_plot("mu_Cs_denser2.png")

# Chern numbers vs mu and BZ grid
# mu, Bz, *Cs = load_data("mu_Bz_Cs.dat")
# Cs = np.sum(Cs[:2], axis=0)
# Z, e = prep_image_data(mu, Bz, Cs)

# # padding do rysowanie dna
# # Z=np.pad(Z,((150, 75),(100,100)), mode='constant', constant_values=np.nan)
# # e=[-1.0, 1.0, -0.5, 0.5]

# plt.imshow(Z, extent=e, origin="lower", aspect="auto", cmap="RdBu", vmin=-2, vmax=2)
# plt.colorbar()
# set_labels(r"$\mu$ (meV)", r"$B_z$ (T)", r"Chern Number for BZ grid $2000 \times 2000$")
# save_plot("mu_Bz_Cs.png")

# # Berry Curvature
# for i, filename in enumerate(["BC.dat", "BC_num.dat", "BC_num2.dat"]):
#     kx, ky, *BCs = load_data(filename)
#     fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(5, 5))


#     for iBC, BC in enumerate(BCs):
#         # print(iBC, np.sum(BC))
#         BC = np.sign(BC) * np.log10((np.abs(BC) + 1.0))
#         Z, e = prep_image_data(kx, ky, BC)
#         cmax = np.abs(BC).max()
#         mapped = axes.flat[iBC].imshow(
#             Z, extent=e, origin="lower", aspect="equal", cmap="RdBu", vmin=-cmax, vmax=cmax
#         )
#         axes.flat[iBC].set_xlabel(r"$k_x$")
#         axes.flat[iBC].set_ylabel(r"$k_y$")
#         axes.flat[iBC].label_outer()
#         axes.flat[iBC].set_title(f"Band {iBC}")
#         axes.flat[iBC].set_box_aspect(1)
#         fig.colorbar(mapped, ax=axes.flat[iBC])

#     fig.suptitle(r"SymLog10 of Berry Curvature for $\mu = -2.89$ meV")
#     save_plot(f"{filename[:-4]}.png")

# # Gap as a function of momentum
# kx, ky, *gaps = load_data("gap.dat")
# kxu, kyu = np.unique(kx), np.unique(ky)
# gap = gaps[0]

# kxe, kye, *energies = load_data("BS.dat")
# kxe, kye = np.meshgrid(np.unique(kx), np.unique(ky), indexing="ij")

# Z, e = prep_image_data(kx, ky, gap)
# plt.imshow(Z, extent=e, origin="lower", aspect="equal", cmap='jet' , vmin=0)
# plt.colorbar()

# set_labels(r"$k_x$", r"$k_y$", r"Gap (meV)")
# save_plot("gap.png")

# # for gap interpolation
# import scipy as scp
# interp = scp.interpolate.RegularGridInterpolator((kxu, kyu), np.reshape(gap, (len(kxu), len(kyu))))

# #for plotting gap along contour
# import matplotlib.cm as cm
# import matplotlib.colors as mcolors

# xyzs = []

# for i,e in enumerate(energies):
#     co=plt.contour(kxe, kye, e.reshape(len(kxe), len(kye)), levels=[0], colors="black", linewidths=0.5)
    
#     # retrive gap values along the contour
#     for poly in co.allsegs[0]:
#         if len(poly) < 3:
#             continue
        
#         x = poly[:,0]
#         y = poly[:,1]
#         z=interp(poly)
        
#         xyzs.append((x,y,z))
        
# plt.close("all")

# all_z = np.concatenate([z for _,_,z in xyzs])
# norm = mcolors.Normalize(vmin=all_z.min(), vmax=all_z.max())
# cmap = cm.jet  # Colormap

# for x,y,z in xyzs:
#     colors = cmap(norm(z))  # Map z-values to colors
#     for i in range(len(x) - 1):
#         plt.plot(x[i:i+2], y[i:i+2], color=colors[i], linewidth=2)

# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# plt.colorbar(sm, ax=plt.gca())
        
# set_labels(r"$k_x$", r"$k_y$", r"Gap (meV)")

# save_plot("gap_contour.png")

# for i,(x,y,z) in enumerate(xyzs):
#     theta = np.arctan2(y,x)
#     idx = np.argsort(theta)
#     plt.plot(theta[idx], z[idx], linewidth=1)
#     set_labels(r"$\theta$", r"Gap (meV)")
#     plt.xticks(np.linspace(-np.pi, np.pi, 5), [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
#     save_plot(f"gap_contour_{i}.png")

# for ax in axes.flat:
    # ax.set_xlabel(r"$\theta$")
    # ax.set_ylabel(r"Gap (meV)")
    # ax.set_yscale("log")
    # ax.set_xticks(np.linspace(-np.pi, np.pi, 5))
    # ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
    # ax.set_ylim(1e-1, 1e-0)
    # ax.set_ylim(1e-1, 1e-0)
#     ax.label_outer()
# fig.tight_layout()
# fig.savefig("plot/gap_contour.png", )
    
# ## FS contour
# n_contours=1
# for i in range(n_contours):
#     kx,ky, *_ = load_data(f"gap_FS{i}.dat")
#     plt.plot(kx, ky, label=f"{i}")

# set_labels(r"$k_x$", r"$k_y$", r"FS")
# # plt.xlim(-0.3, 0.3)
# # plt.ylim(-0.3, 0.3)
# plt.gca().set_box_aspect(1)
# # plt.legend(title="contour idx")
# save_plot("FS.png")

# ## Gap along FS contour
# kx_ky_gap_lst = []
# for i in range(n_contours):
#     kx,ky, *gaps = load_data(f"gap_FS{i}.dat")
#     kx_ky_gap_lst.append((kx, ky, gaps[0]))
#     theta=np.arctan2(ky[1:-1], kx[1:-1])
#     gap = gaps[0][1:-1]
#     plt.plot(theta, gap, label=f"{i}")

# # plt.yscale("log")
# plt.xlim(-np.pi, np.pi)
# plt.xticks(np.linspace(-np.pi, np.pi, 5), [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
# # plt.legend(title="contour idx")
# set_labels(r"$\theta$", r"$\Delta$", r"Gap along FS contours")
# save_plot("gap_FS.png")

# # #for plotting gap along contour
# import matplotlib.cm as cm
# import matplotlib.colors as mcolors

# all_g = np.concatenate([g for _,_,g in kx_ky_gap_lst])
# norm = mcolors.Normalize(vmin=all_g.min(), vmax=all_g.max())
# cmap = cm.jet  # Colormap

# for x,y,z in kx_ky_gap_lst:
#     colors = cmap(norm(z))  # Map z-values to colors
#     for i in range(len(x) - 1):
#         plt.plot(x[i:i+2], y[i:i+2], color=colors[i])

# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# plt.colorbar(sm, ax=plt.gca())
    
# set_labels(r"$k_x$", r"$k_y$", r"Gap (meV)")

# save_plot("gap_colored_FS.png")