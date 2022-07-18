
function time_step_iteration( curr_solution, model_properties; verbosity=false )

    curr_time, Un, Udotn, Udotdotn, system_matrix, system_rhs = unwrap( curr_solution )

    # start iteration
    converged_boolean = false
    iters = 0

    # valid for the static case only !!!
    Unp1k  = Un 
    #
    Udotnp1    = Udotn
    Udotdotnp1 = Udotdotn

    while converged_boolean == false
        
        # add to iterations counter
        iters = iters + 1

        # compute deltaU
        deltau_red = compute_DeltaU( system_matrix, system_rhs )

        # update Uk
        Unp1k[ neum_dofs] = Unp1k[ neum_dofs ] + deltau_red

        # compute system
        system_matrix, system_RHS = assemble_system( Unp1k, model_properties )

        # convergence test
        booleanConverged = convergenceTest( modelProperties.analysisSettings, [], FextG(BCsData.neumDofs), deltaured, Utp1k(BCsData.neumDofs), iters, [], systemDeltauRHS ) ;

        # print

    end
    next_solution = curr_solution
    next_solution.time = curr_solution.time + model_properties.analysis_settings.delta_time

    return next_solution
end




function compute_DeltaU( matrix, rhs )
    DeltaU = matrix \ rhs
    return DeltaU
end



function convergence_test( analysis_settings, redFint, redFext, redDeltaU, redUk, dispIter, redFinet, systemDeltauRHS )
  
    stop_tol_disps = analysis_settings.stop_tol_disps
    stop_tol_force = analysis_settings.stop_tol_force
    stop_tol_iters = analysis_settings.stop_tol_iters   
  
    # normaUk       = norm( redUk )               
    # normadeltau   = norm( redDeltaU         )   
  
    # if strcmp( analysisSettings.methodName, "arcLength")
    #   systemDeltauRHS = systemDeltauRHS(:,1);
    # end
  
    # deltaErrLoad  = norm( systemDeltauRHS )     ;
    # normFext      = norm( redFext         )     ;
  
    # logicDispStop = ( normadeltau  < ( normaUk  * stopTolDeltau ) )  ;
    # logicForcStop = ( deltaErrLoad < ( (normFext+(normFext < stopTolForces)) * stopTolForces ) )  * ( deltaErrLoad > 0 ) ;
  
    # if logicForcStop
    #   stopCritPar = 1 ;      booleanConverged = 1 ;
    # elseif logicDispStop
    #   stopCritPar = 2 ;      booleanConverged = 1 ;
    # elseif ( dispIter >= stopTolIts )
    #   stopCritPar = 3 ;      booleanConverged = 1 ;
    # else
    #   stopCritPar = [];      booleanConverged = 0 ;
    # end

    return booleanConverged, stopCritPar, deltaErrLoad
end
    
