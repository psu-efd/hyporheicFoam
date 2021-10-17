 
Info<< "Reading surface field surYi\n" << endl;
  while (pimple.loop())
   {
     if(CSurfSwitch)
      { 
         surreaction->correct();
         volScalarField surYi //surface scalar 
         (
             IOobject
             (
                forAll(surY, i)
                {
                    if (surcomposition.active(i))
                         {
                             surYi =surY[i];
                             runTime.timeName(),
                             meshSurface,
                             IOobject::MUST_READ,
                             IOobject::AUTO_WRITE
                         }
                } 
             ),
             meshSurface
          ),
      }
    };
