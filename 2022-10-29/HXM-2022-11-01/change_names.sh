# After command 'gmshIO -gmsh_file ...', we got some vtu/vtp files.
# Change their names by refering to '../PERIGEE/examples/tet4_fsi/preprocess_fsi.cpp'.
# 'f' means 'fluid' ------> lumen
# 's' means 'solid' ------> tissue
# 'f_s' ------> interior
# KEEP 'whole_vol.vtu'.

mv fluid.vtu lumen_vol.vtu
mv fwall_fluid.vtp lumen_wall_vol.vtp
mv fbot_fluid.vtp lumen_inlet_vol_000.vtp
mv ftop_fluid.vtp lumen_outlet_vol_000.vtp

mv fwall_solid.vtp tissue_interior_wall_vol.vtp

mv solid.vtu tissue_vol.vtu
mv swall_solid.vtp tissue_wall_vol.vtp
mv sbot_solid.vtp tissue_inlet_vol_000.vtp
mv stop_solid.vtp tissue_outlet_vol_000.vtp


