/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createTimeControls.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    Foam::clock clockTime;

    Random rndGen(clockTime.getTime());
    // Gen random seed
    // rndGen.seed();
    const volVectorField& C = mesh.C();

    #include "createFields.H"

    dimensionedScalar velT = fvc::domainIntegrate(DK * T*(1-T)) ;
    dimensionedScalar vel = fvc::domainIntegrate(DK * T*(1-T)) ;
    dimensionedScalar velRes (fvc::domainIntegrate(Db*xi - DK * T*xi*xi));
    // scalar TInitialResidual(1.0);
    // scalar xiInitialResidual(1.0);
    label icount(0);
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    // // Update velocity

    dimensionedScalar tvel(0.0);
    // dimensionedScalar dsigma = Foam::sqrt(2.0e0 * dtheta);

    Info << "Theta = " << dtheta.value() << " Sigma = " << dsigma.value()<<endl;

    
    dimensionedScalar dx = Foam::cmptMag(C[1][0] - C[0][0]);
    scalar rcorr(Foam::exp(-dx.value()*dtheta.value()));
    scalar rroot(Foam::sqrt(1.0-rcorr*rcorr));
    scalar dsigma_t = dsigma.value() / Foam::sqrt(2.0 * dtheta.value());
    Info << "dx = " << dx.value() << " corr = " << rcorr << " sqrt(r) = "<< rroot<< nl << endl;
    List<scalar> dF(C.size());
    List<scalar> dr(C.size());

    if (!useOldField) {
      Info << "Generate a new stochastic field" << nl << endl;
    forAll(U, i){
        dimensionedScalar dW = rndGen.scalarNormal();
        dimensionedScalar dWi = dW * dsigma_t;
        // dW = dW * dsigma * Foam::sqrt(dx.value());
        dW *= Foam::sqrt(dx.value());
        dr[i] = dW.value();
        if (i == 0){
 //           U[i][0] = dW.value() ;
            dF[i] = dWi.value();
        } else {
            //U[i][0] = U[i-1][0] + dtheta.value()*(dmean.value() -U[i-1][0] ) * dx.value() + dW.value();
            // U[i][0] = (1.0-dtheta.value()*dx.value()) * U[i-1][0] + dW.value();
            dF[i] = rcorr * dF[i-1] + rroot * dWi.value();
        }
        velInit[i] = dF[i];
        velInit[i] = dW.value();

        
//        Info << "Mesh " << C[i][0] << " " << U[i][0] << endl;
      }
      // calc mean value of velInit and modify it to mean zero
      scalar velm(0.0);
      forAll(velInit, i){
          velm += velInit[i];
      }
      velm /= velInit.size();
      forAll(velInit, i){
          velInit[i] -= velm;
      }
    } else {

      Info << "Use existing stochastic field" << nl << endl;

    }
    // volScalarField kExpWn(Foam::exp(-Foam::sqrt(barVel.value())*velInit));
    // volScalarField kExpWp(Foam::exp(Foam::sqrt(barVel.value())*velInit)*DT);
    volScalarField kExpWn(Foam::exp(-Foam::sqrt(barVel.value())*velInit));
    volScalarField kExpWp(Foam::exp(Foam::sqrt(barVel.value())*velInit)*DT);
    forAll(velInit, i){
        // scalar tt = Foam::sqrt(barVel.value())*velInit[i];
        scalar tt = barVel.value() * velInit[i];
        kExpWn[i] = Foam::exp(-tt);
        kExpWp[i] = Foam::exp(-tt)*DT.value();
    }

    // Scale the values by barVel
    //
    // scalar vShift = barVel.value()*2.0*dsigma_t;
    scalar velm(0.0);
    scalar vels(0.0);
    forAll(U,i){
        U[i][0] = velInit[i] * barVel.value();
        U[i][0] -= vShift.value();
//         U[i][0] = 2.0;
      //  U[i][0] = 0.0;
      velm += U[i][0];
    }
    velm /= mesh.C().size();
    forAll(U,i) {
      vels += (U[i][0] - velm)*(U[i][0]-velm);
    }
    vels /= mesh.C().size();
    Info << "VelField "<<velm <<" "<< vels << endl;
    std::ofstream file;
    // file.open("res.dat",std::ofstream::out| std::ofstream::app);
    file.open("res.dat",std::ofstream::out);
    forAll(dF, i) {
      file << i<<","<< C[i][0] << ',' << velInit[i]<<","<<U[i][0] << ',' << dr[i]<<std::endl;
    }
    file.close();

    #include "setInitialDeltaT.H"

    Info << "Start timesteps" << endl;

    while (simple.loop(runTime))
    {

        // forAll(U, i){
        //     scalar dW = rndGen.scalarNormal();
        //     velInit[i] = dW;
        //     U[i][0] = velInit[i] * barVel.value();
        // }

        #include "readTimeControls.H"

        #include "CourantNo.H"

        #include "setDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

//             // Update velocity
//             vel = fvc::domainIntegrate(DK * T*(1-T)) ;
//             velT = fvc::domainIntegrate(DK * T*T*(1-T))       ; 
//             forAll(U,i){
//                 U[i][0] = vel.value();
//             }
        // Update phi
        //
        phi = fvc::flux(U);

        volDivPhi = fvc::div(phi);
        volDivU = fvc::div(U);
        
        icount = 0;
        bool cond(true);
        // while (simple.correctNonOrthogonal())
        // while (((xiInitialResidual > 1.0e-6 )||(TInitialResidual > 1.0e-6 ) ) && (icount<5))
        while (cond)
        {

            icount += 1;
            // fvScalarMatrix TEqn
            // (
            //     fvm::ddt(T)
            //   - fvm::div(phi, T)
            //   - fvm::laplacian(DT, T)
            //   + fvm::Sp(2*DK*T + fvc::div(U) - DK, T)
            //   ==
            //     // fvm::Sp(DK - DK*T - fvc::div(U), T)
            //     // fvm::Sp(DK - 2*DK*T - fvc::div(U), T)
            // //   + fvc::Su(DK*T*T,T)
            //     fvc::Su(DK*T*T,T)
            // );

            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              - kExpWp * fvm::laplacian(kExpWn, T)
              + fvm::Sp(2*DK*T - DK, T)
              ==
                fvc::Su(DK*T*T,T)
            );

            TEqn.relax();
            fvOptions.constrain(TEqn);
            // TEqn.solve();
            scalar TInitialResidual = TEqn.solve().initialResidual();
            fvOptions.correct(T);

            /*fvScalarMatrix xiEqn
            (
                fvm::ddt(xi)
              + fvm::div(phi, xi)
            //   - T * fvc::div(phi)
              - fvm::laplacian(Dxi, xi)
            //   - fvm::SuSp(fvc::div(phi)+DK*T*xi-Db, xi)
              - fvm::Sp(fvc::div(phi)+DK*T*xi-Db, xi)
             ==
                fvOptions(xi)
            );

            xiEqn.relax();
            fvOptions.constrain(xiEqn);
            scalar xiInitialResidual = xiEqn.solve().initialResidual();
            fvOptions.correct(xi);
            */
            // TEqn.solve();
            // xiEqn.solve();
            scalar xiInitialResidual(0.0);
            // Info<< "Residuals = " << TEqn.solve().max().initialResidual() << " " << xiEqn.solve().max().initialResidual()  << endl;
            Info<< "Residuals = " << TInitialResidual << " " << xiInitialResidual << " Cnt = "<< icount  << endl;

            if (TInitialResidual<1.0e-5 && xiInitialResidual < 1.0e-5) {
                cond = false ;
            }
            if (icount>= 25) {
                cond = false;
            }
        }

        // dimensionedScalar vel(fvc::domainIntegrate(T*T*(1-T))) ;
        vel = fvc::domainIntegrate(DK * T*(1-T)) ;
        velT = fvc::domainIntegrate(DK * T*T*(1-T)) ;
        velRes = fvc::domainIntegrate(Db*xi - DK * T*xi*xi) ;
        scalar flamePos(0.0);
        scalar meanVel(0.0);
        forAll(mesh.C(),i) {
            meanVel += U[i][0];
        }
        forAll(mesh.C(),i) {
            if (T[i] >= 0.50e0) {
                flamePos = C[i].component(0);
                break;
            }
        }
        meanVel /= mesh.C().size();
        Info<< nl << "Velocity = " << vel.value()<< " " << velRes.value() << " " << tvel.value() << " " << flamePos << " " << Foam::gMax(T.internalField()) << " " << Foam::gMax(xi.internalField())<< " " << meanVel << nl << endl;
        // Info<< "phi = " << phi.value()<< " " << velRes.value() << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
