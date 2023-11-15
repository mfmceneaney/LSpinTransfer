import numpy as np

def compute_systematics(results,bin_migration_mat=None,bin_migration_order=1,systematic_scales_mat=None):
    """
    :description:

    :param: result
    :param: bin_migration_mat
    :param: systematic_scales_mat

    :returns: systematics
    """
    results_order = len(results.shape) #NOTE: THIS IS JUST WHAT ORDER BINNING YOU ARE DEALING WITH.
    systematics = np.zeros(results.shape)

    # Compute and add bin migration systematics
    if bin_migration_mat is not None: #NOTE: ASSUME BINNING IS 1D HERE.
        if results_order==1:
            diags = np.ones(np.diag(bin_migration_mat).shape)
            neighbors = np.sum(
                [
                    np.sum((np.diag(diags[i:],i),np.diag(diags[i:],-i)),axis=0)
                    for i in range(1,bin_migration_order+1)
                ],
                axis=0)
            systematics += np.diag( # DeltaA = Diag(n.(a.f-f.a)) = Diag(n.[a,f])
                np.matmul(
                    neighbors,
                    np.matmul(np.diag(results),bin_migration_mat) - np.matmul(bin_migration_mat,np.diag(results))
                )
            )
        else:
            print("WARNING: BIN MIGRATION NOT IMPLEMENTED FOR 2+D BINNING CASES.")

    # Apply multiplicative scale systematics, note that these should already be summed over all sources of systematic error
    if systematic_scales_mat is not None:
        systematics += np.multiply(results,systematic_scales_mat)

    return systematics

results = np.array([0.1, 0.2, 0.0, 0.0, 0.0],dtype=float)
bin_migration_mat = np.array([ #NOTE: (i,j) -> (from_bin,to_bin) so if (i=0,j=1)=0.1 AND results[i=0]=0.1 AND results[i=1]=0.2 AND all else 0
                                                                        # THEN systematics[i=1]= -f[i=0,j=1]*results[i=1] + f[i=0,j=1]*results[i=0]
                                                                        #                      = -0.1       * 0.2         + 0.1       * 0.1
                                                                        #                      = -0.01
    [0.0, 0.1, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0],
],dtype=float)

systematic_scales_mat = None #np.array([0.0, 0.0, 0.1, 0.0, 0.0],dtype=float)

systematics = compute_systematics(
    results,
    bin_migration_mat=bin_migration_mat,
    bin_migration_order=1,
    systematic_scales_mat=systematic_scales_mat
    )
print("DEBUGGING: results = ",results)
print("DEBUGGING: systematics = ",systematics)
