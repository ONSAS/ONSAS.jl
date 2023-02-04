
"""
Function that computes the matrix and Right-hand-side of the system of equations.
"""
function assemble_system(Model_properties::ModelProperties, U::Vector, neum_dofs)


    Algorithm = Model_properties.Algorithm
    BCsData = Model_properties.Boundary_Conditions

    # computes internal forces & matrices
    fs, mats = assembler(Model_properties.Materials, Model_properties.Geometries, Model_properties.Mesh, Sol)

    # assemble user external forces
    Fint = fs[1]
    KT = mats[1]

    # if cmp( analysisSettings.method, "new"' ) || strcmp( analysisSettings.methodName, 'alphaHHT' )
    #     dampingMat = mats{2} ;
    #     massMat    = mats{3} ;

    # if cmp(analysisSettings.method, "newton_raphson") == 0 || cmp(analysisSettings.method, "arc_length") == 0

    #     systemΔuMatrix = KT[neumdofs, neumdofs]

    #     [FextG, nexTimeLoadFactors] = computeFext(BCsData.factorLoadsFextCell, BCsData.loadFactorsFuncCell, modelProperties.analysisSettings, nextTime, length(Fint), BCsData.userLoadsFilename, [])

    #     rhs = -(Fint[BCsData.neumDofs] - FextG[BCsData.neumDofs] - Faero(BCsData.neumDofs))

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

    # end

    systemΔuMatrix, rhs = systemMatrices(Algorithm::AbstractAlgorithm, KT::Matrix, Fint::Vector, BCsData, neum_dofs)

    return systemΔuMatrix, rhs
end


function systemMatrices(Algorithm::NewtonRaphson, KT::MAtrix, Fint::Vector, BCsData, neum_dofs)

    # Tangent stiffness matrix
    systemΔuMatrix = KT[neum_dofs, neum_dofs]
    # External force vector
    FextG = zeros(length(Fint))
    FextG = compute_fext!(Algorithm::AbstractAlgorithm, factorLoadsFextCell, loadFactorsFuncCell, userLoadsFilename, time, Fext, [])
    # Residual
    # rhs = -(Fint[BCsData.neumDofs] - FextG[BCsData.neumDofs] - Faero(BCsData.neumDofs))
    rhs = -(Fint[neum_dofs] - FextG[neum_dofs])

    return systemΔuMatrix, rhs
end