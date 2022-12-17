
"""
Function that computes the matrix and Right-hand-side of the linear system of equations.
"""
function assemble_system(Sol::ModelSolution, Model_properties::ModelProperties)

    Analysis_settings = Model_properties.analysis_settings

    # computes internal forces
    fs, mats = assembler(Model_properties.materials, Model_properties.geometries, Model_properties.mesh, Sol)

    # assemble user external forces
    Fint = fs[1]
    Kint = mats[1]

    # if cmp( analysisSettings.method, "new"' ) || strcmp( analysisSettings.methodName, 'alphaHHT' )
    #     dampingMat = mats{2} ;
    #     massMat    = mats{3} ;

    if cmp(analysisSettings.method, "newton_raphson") == 0 || cmp(analysisSettings.method, "arc_length") == 0

        matrix = Kint[neumdofs, neumdofs]

        [FextG, nexTimeLoadFactors] = computeFext(BCsData.factorLoadsFextCell, BCsData.loadFactorsFuncCell, modelProperties.analysisSettings, nextTime, length(Fint), BCsData.userLoadsFilename, [])

        rhs = -(Fint[BCsData.neumDofs] - FextG[BCsData.neumDofs] - Faero(BCsData.neumDofs))

        # elseif strcmp( analysisSettings.methodName, 'newmark' )

        #   alphaNM = analysisSettings.alphaNM ;
        #   deltaNM = analysisSettings.deltaNM ;
        #   deltaT  = analysisSettings.deltaT  ;

        #   systemDeltauMatrix =                                   KT(         neumdofs, neumdofs ) ...
        #                        + 1/( alphaNM * deltaT^2)       * massMat(    neumdofs, neumdofs ) ...
        #                        + deltaNM / ( alphaNM * deltaT) * dampingMat( neumdofs, neumdofs )  ;

        # elseif strcmp( analysisSettings.methodName, 'alphaHHT' )

        #   alphaHHT = analysisSettings.alphaHHT ;
        #   deltaT   = analysisSettings.deltaT  ;

        #   deltaNM = (1 - 2 * alphaHHT ) / 2 ;
        #   alphaNM = (1 - alphaHHT ^ 2 ) / 4 ;

        #   systemDeltauMatrix = (1 + alphaHHT )                                 * KT         ( neumdofs, neumdofs ) ...
        #                      + (1 + alphaHHT ) * deltaNM / ( alphaNM*deltaT  ) * dampingMat ( neumdofs, neumdofs )  ...
        #                      +                         1 / ( alphaNM*deltaT^2) * massMat    ( neumdofs, neumdofs ) ;

    end

    return matrix, rhs
end