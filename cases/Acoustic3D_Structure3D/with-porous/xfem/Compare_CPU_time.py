import scipy

def h2s(h,m,s):
    return h*3600+m*60+s


# On one proc only.
# 40 freq. step

##Number of nodes: 29177
##Number of elements in air: 148683
##Number of elements in porous: 4466
##Number of nodes in air: 27068
##Number of nodes in porous: 1369
##Number of nodes at interface: 555
##nnodes for structure= 3275
##nelem for structure= 6335

# NO REDUCTION / NO POROUS: 40 steps:
# Proc.  0  / time at the beginning of the FRF: Mon Apr 24 17:05:59 2017
# Proc.  0  / time at the end of the FRF: Mon Apr 24 17:08:29 2017
online_ref_no_porous=h2s(17,8,29)-h2s(17,5,59)

##########################################################
# XFEM, no reduction : Main_xfem_fluid_no_porous_flex_struc_no_reduction_8.py
# time at the beginning of the computation: Mon Apr 24 16:53:55 2017
# Proc.  0  / time at the beginning of the FRF: Mon Apr 24 16:54:12 2017
# Proc.  0  / time at the end of the FRF: Mon Apr 24 16:57:12 2017
# nb of total dofs: 50571

offline_ref=h2s(16,54,12)-h2s(16,53,55)
online_ref=h2s(16,57,12)-h2s(16,54,12)

##########################################################
# XFEM, reduction : Main_xfem_fluid_porous_flex_struc_CB_reduction_8.py
# 140F+100S modes
# time at the beginning of the computation: Mon Apr 24 16:57:12 2017
# Proc.  0  / time at the beginning of the FRF: Mon Apr 24 16:59:30 2017
# Proc.  0  / time at the end of the FRF: Mon Apr 24 17:00:45 2017
# Last eigenfrequency in fluid basis:  234.696045066
# Last eigenfrequency in structure basis:  310.121797212
# nb of total dofs: 5105

offline_reduced=h2s(16,59,30)-h2s(16,57,12)
online_reduced=h2s(17,00,45)-h2s(16,59,30)


print(' CASE           /  NON-REDUCED  /  REDUCED / NON-REDUCED, NO POROUS ')
print('OFFLINE             ',offline_ref    ,'              ',offline_reduced)
print('ONLINE, 40 steps    ',online_ref     ,'              ',online_reduced,'           ',online_ref_no_porous)
print('ONLINE, 1 step      ',online_ref/40.0,'         ',online_reduced/40.0,'           ',online_ref_no_porous/40.0)

#################################
## FINE MESH
#################################
online_ref_no_porous=h2s(19,9,6)-h2s(19,1,2)

offline_ref=h2s(18,42,4)-h2s(18,41,44)
online_ref=h2s(18,51,56)-h2s(18,42,4)

offline_reduced=h2s(18,55,59)-h2s(18,52,5)
online_reduced=h2s(19,0,37)-h2s(18,55,59)

print('')
print('')
print('     FINE MESH')
nbsteps= 20.0*80.0
print(' CASE                /  NON-REDUCED  /  REDUCED / NON-REDUCED, NO POROUS ')
print('OFFLINE                  ',offline_ref    ,'              ',offline_reduced)
print('ONLINE, ',nbsteps,' steps    ',online_ref     ,'              ',online_reduced,'           ',online_ref_no_porous)
print('ONLINE, 1 step           ',online_ref/nbsteps,'         ',online_reduced/nbsteps,'           ',online_ref_no_porous/nbsteps)

#################################
## CONVERGED MESH
#################################

##One proc. / 40 steps
##
##Number of nodes: 21688
##Number of elements in air: 108220
##Number of elements in porous: 3559
##Number of nodes in air: 19972
##Number of nodes in porous: 1112
##Number of nodes at interface: 453
##nnodes for structure= 3275
##nelem for structure= 6335
##nb_mode_F = 200
##nb_mode_S = 70
##Last eigenfrequency in fluid basis:  273.656308676
##Last eigenfrequency in structure basis:  216.159494654

offline_reduced_0=h2s(13,21,51)-h2s(13,20,59)
offline_reduced_1=h2s(13,25,9)-h2s(13,24,32)
online_reduced=h2s(13,25,52)-h2s(13,25,9)
#nb of total dofs:  4078

offline_ref=h2s(13,27,57)-h2s(13,27,45)
online_ref=h2s(13,29,40)-h2s(13,27,57)
#nb of total dofs:  42520

print('\n')
print('\n')
print('--------------------------')
print('One proc. / 40 steps')
print('--- NO REDUCTION ---')
print('Offline : ',offline_ref)
print('Online 40 steps : ',online_ref)
print('Online 1 step : ',online_ref/40.0)
print('Total, 400 steps : ',offline_ref+400.0*online_ref/40.0)
print('\n')
print('--- REDUCTION ---')
print('Offline 0 : ',offline_reduced_0)
print('Offline 1 : ',offline_reduced_1)
print('Online 40 steps : ',online_reduced)
print('Online 1 step : ',online_reduced/40.0)
print('Total, 400 steps : ',offline_reduced_1+400.0*online_reduced/40.0)

print('gain factor, online 1 step : ',online_ref/online_reduced)

print('gain factor, 400 steps : ',(offline_ref+400.0*online_ref/40.0)/(offline_reduced_1+400.0*online_reduced/40.0))


##20 proc. / 80 steps
print('\n')
print('\n')
print('--------------------------')
print('20 proc. , 80 steps/proc ==> 1600 steps')
print('--- NO REDUCTION ---')
offline_reduced_0=h2s(14,10,47)-h2s(14,9,49)
offline_reduced_1=h2s(13,52,52)-h2s(13,52,14)
online_reduced=h2s(13,54,46)-h2s(13,52,52)
#nb of total dofs:  4078

offline_ref=h2s(13,58,46)-h2s(13,58,31)
online_ref=h2s(14,3,26)-h2s(13,58,46)
print('--- NO REDUCTION ---')
print('Offline : ',offline_ref)
print('Online 20*80 steps : ',online_ref)
print('Online 1 step, 1 proc : ',online_ref/80.0)
print('Online 1 step, 20 proc : ',online_ref/1600.0)
print('\n')
print('--- REDUCTION ---')
print('Offline 0 (mumps, 20 proc) : ',offline_reduced_0)
print('Offline 1 : ',offline_reduced_1)
print('Online 20*80 steps : ',online_reduced)
print('Online 1 step, 1 proc : ',online_reduced/80.0)
print('Online 1 step, 20 proc : ',online_reduced/1600.0)
