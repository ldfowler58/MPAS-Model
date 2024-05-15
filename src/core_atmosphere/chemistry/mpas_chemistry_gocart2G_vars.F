!=================================================================================================================
 module mpas_chemistry_gocart2G_vars
!=================================================================================================================
 use CA2G_bc_GridCompMod
 use CA2G_br_GridCompMod
 use CA2G_oc_GridCompMod
 use DU2G_GridCompMod
 use NI2G_GridCompMod
 use SS2G_GridCompMod
 use SU2G_GridCompMod

 use CA2G_bc_StateSpecs,only: CA2G_bc_State
 use CA2G_br_StateSpecs,only: CA2G_br_State
 use CA2G_oc_StateSpecs,only: CA2G_oc_State
 use DU2G_StateSpecs,only: DU2G_State
 use NI2G_StateSpecs,only: NI2G_State
 use SU2G_StateSpecs,only: SU2G_State

 use mpas_chemistry_gocart2G_interface,only: atm_gocart2G


 type(CA2G_bc_GridComp):: CA2G_bc_params
 type(CA2G_br_GridComp):: CA2G_br_params
 type(CA2G_oc_GridComp):: CA2G_oc_params
 type(DU2G_GridComp):: DU2G_params
 type(NI2G_GridComp):: NI2G_params
 type(SS2G_GridComp):: SS2G_params
 type(SU2G_GridComp):: SU2G_params

 type(CA2G_bc_State):: CA2G_bc
 type(CA2G_br_State):: CA2G_br
 type(CA2G_oc_State):: CA2G_oc
 type(DU2G_State):: DU2G
 type(NI2G_State):: NI2G
 type(SU2G_State):: SU2G

 type(atm_gocart2G):: mpas_gocart2G

!=================================================================================================================
 end module mpas_chemistry_gocart2G_vars
!=================================================================================================================

