 #!/usr/bin/env python
# encoding: utf-8

from matplotlib import rc
rc("font", family="serif", size=8)
rc("text", usetex=True)

# find this!
import daft

# Colors for different distributions! see Wiecki paper
gamma_color = {"ec": "#46a546"}
normal_color = {"ec": "#f89406"}
halfnormal_color = {"ec": "#0000cd"}
invlogit_color = {"ec": "#d02090"}

pgm = daft.PGM([6,6])
pgm.add_node(daft.Node("x", r"$x_{(p,s,t)}$", 3, 3.5, observed=True))
 # plates - xstart, ystart, xsize, ysize
pgm.add_plate(daft.Plate([2.45, 3, 1.1, 1], label=r"$t = 1, ..., $trials"))

pgm.add_node(daft.Node("z", r"$z_{(p,s)}$", 2, 3.2, plot_params=invlogit_color))
pgm.add_edge("z", "x")
pgm.add_node(daft.Node("mz", r"$\mu_{z(p)}$", 1, 3.2, plot_params=normal_color))
pgm.add_edge("mz", "z")
pgm.add_node(daft.Node("sz", r"$\sigma_{z}$", 1, 2.5, plot_params=halfnormal_color))
pgm.add_edge("sz", "z")

pgm.add_node(daft.Node("dc", r"$dc_{(p,s)}$", 2, 3.8, plot_params=normal_color))
pgm.add_edge("dc", "x")
pgm.add_node(daft.Node("mdc", r"$\mu_{dc(p)}$", 1, 3.8, plot_params=normal_color))
pgm.add_edge("mdc", "dc")
pgm.add_node(daft.Node("sdc", r"$\sigma_{dc}$", 1, 4.5, plot_params=halfnormal_color))
pgm.add_edge("sdc", "dc")

# plates - xstart, ystart, xsize, ysize
pgm.add_plate(daft.Plate([0.3, 2.8, 3.3, 1.4], label=r"$p = -1, 1$ previous choices"))
     
# things on the right hand side

pgm.add_node(daft.Node("v", r"$v_{(s)}$", 4, 4.1, plot_params=normal_color))
pgm.add_edge("v", "x")
pgm.add_node(daft.Node("mv", r"$\mu_{v}$", 5, 5, plot_params=normal_color))
pgm.add_edge("mv", "v")
pgm.add_node(daft.Node("sv", r"$\sigma_{v}$", 5, 4.4, plot_params=halfnormal_color))
pgm.add_edge("sv", "v")

pgm.add_node(daft.Node("a", r"$a_{(s)}$", 4, 3.5, plot_params=gamma_color))
pgm.add_edge("a", "x")
pgm.add_node(daft.Node("ma", r"$\mu_{a}$", 5, 3.8, plot_params=gamma_color))
pgm.add_edge("ma", "a")
pgm.add_node(daft.Node("sa", r"$\sigma_{a}$", 5, 3.2, plot_params=halfnormal_color))
pgm.add_edge("sa", "a")

pgm.add_node(daft.Node("t", r"$t_{(s)}$", 4, 2.9, plot_params=normal_color))
pgm.add_edge("t", "x")
pgm.add_node(daft.Node("mt", r"$\mu_{t}$", 5, 2.6, plot_params=gamma_color))
pgm.add_edge("mt", "t")
pgm.add_node(daft.Node("st", r"$\sigma_{t}$", 5, 2, plot_params=halfnormal_color))
pgm.add_edge("st", "t")

pgm.add_plate(daft.Plate([1.5, 2.5, 3, 1.9], label=r"$s = 1, ..., $subjects"))

pgm.add_node(daft.Node("svv", r"$sv$", 3, 4.8, plot_params=halfnormal_color))
pgm.add_edge("svv", "x")

# Render and save.
pgm.render()
pgm.figure.savefig("HDDM.pdf")
pgm.figure.savefig("HDDM.png", dpi=150)