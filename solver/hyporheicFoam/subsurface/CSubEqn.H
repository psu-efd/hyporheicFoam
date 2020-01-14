while (simpleSubsurface.correctNonOrthogonal())
{
    fvScalarMatrix CSubEqn
        (
              fvm::ddt(CSub)
            + fvm::div(phiSub, CSub)
            - fvm::laplacian(DSub_Coeff, CSub)
            ==
              fvOptions(CSub)
        );

    CSubEqn.relax();

    fvOptions.constrain(CSubEqn);

    CSubEqn.solve();

    fvOptions.correct(CSub);

    //clip the negative value of concentration
    CSub.max(0.0);
}
