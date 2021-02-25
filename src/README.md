# Femocs Catalyst

The Catalyst adaptor for [FEMOCS](https://github.com/veskem/femocs) library for ParaView prior to v5.9 in which Catalyst API has been changed.

To compile the library and the demo program run:

```shell
make
```

The following library should be installed or compiled:

- Deal.II 9.2
- ParaView 5.8 with the Catalyst flag
- FEMOCS

To specify paths to compiled libraries, pass it to `make`:

```shell
make FEMOCS_DIR=<femocs_dir> PARAVIEW_DIR=<paraview_dir> DEAL_II_DIR=<deal_ii_dir>
```

To clean the build, run:

```shell
make clean
```