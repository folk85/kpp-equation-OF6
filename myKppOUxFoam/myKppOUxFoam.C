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
#include <random>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createTimeControls.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    Foam::clock clockTime;

    std::random_device rd;

    // Random rndGen(clockTime.getTime());
    Random rndGen(rd());
    // Gen random seed
    // rndGen.seed();
    const volVectorField& C = mesh.C();

    #include "createFields.H"

    Info << "Read the sdeScheme: " << sdeScheme << endl;
    word wSdeType("hello");
    if ( sdeScheme == wSdeType ){
        Info << " It says Hello! " << endl;
    }

    dimensionedScalar velT = fvc::domainIntegrate(DK * T*(1-T)) ;
    dimensionedScalar vel = fvc::domainIntegrate(DK * T*(1-T)) ;
    dimensionedScalar velRes (fvc::domainIntegrate(Db*xi - DK * T*xi*xi));
    dimensionedScalar vel_mean = velT;
    dimensionedScalar vel_meano = velT;
    // scalar TInitialResidual(1.0);
    // scalar xiInitialResidual(1.0);
    label icount(0);
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    // // Update velocity

    dimensionedScalar tvel(0.0);
    // dimensionedScalar dsigma = Foam::sqrt(2.0e0 * dtheta);

    Info << "sdeScheme = " << sdeScheme << endl ;
    Info << "theta = " << dtheta.value() << endl ;
    Info << "tau = " << dtau.value() << endl ;
    Info << "barVel = " << barVel.value() << endl ;
    Info << "sdeScheme = " << sdeScheme << endl ;
    Info << "sigma = " << dsigma.value() << endl ;
    Info << "deps = " << deps << endl ;
    Info << "dnu = " << dnu << endl ;

    // scalar timeo = runTime.time();
    // dimensionedScalar dx = Foam::cmptMag(C[1][0] - C[0][0]);
    scalar dx = Foam::cmptMag(C[1][0] - C[0][0]);
    scalar rcorr(Foam::exp(-dx*dtheta.value()));
    scalar rroot(Foam::sqrt(1.0-rcorr*rcorr));
    scalar dsigma_t = dsigma.value() / Foam::sqrt(2.0 * dtheta.value());
    Info << "dx = " << dx << " corr = " << rcorr << " sqrt(r) = "<< rroot<< nl << endl;
    List<scalar> dF(C.size());
    List<scalar> dFo(C.size());
    List<scalar> dr(C.size());
    label nel = C.size();
    scalar dWs(0.0);
    scalar dWi(0.0);


    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::normal_distribution<double> dnd{0,1};

    if (!useOldField) {
      Info << "Generate a new stochastic field" << nl << endl;

      // Set Stochastic field Orstein-Uhlenbeck along the X-space as an initial conditions
      // It works with the old version of cases. If the variable "sdeScheme" isn't set in transportProperties
      // then the program uses the O-U process. 
      //
      // #include "gen_stoch_field.H"

      // word wSdeOU("Orstein");
      // if (sdeScheme == wSdeOU ) {
      if ((sdeScheme == word("Orstein"))||(sdeScheme == word("OrsteinTime") )) {
        dsigma_t = Foam::sqrt(0.5e0 * dtheta.value());
        forAll(U, i){
          // dimensionedScalar dW = rndGen.scalarNormal();
          dWs = rndGen.scalarNormal();
          // dWs = dnd(gen);
          dWi = dWs * dsigma_t;
          // dW = dW * dsigma * Foam::sqrt(dx);
          // dW /= Foam::sqrt(dx);
          dr[i] = dWs;
          if (i == 0){
  //           U[i][0] = dW.value() ;
              dF[i] = dWi;
          } else {
              //U[i][0] = U[i-1][0] + dtheta.value()*(dmean.value() -U[i-1][0] ) * dx + dW.value();
              // U[i][0] = (1.0-dtheta.value()*dx) * U[i-1][0] + dW.value();
              dF[i] = rcorr * dF[i-1] + rroot * dWi;
          }
          velInit[i] = dF[i];
          // velInit[i] = dW.value();
        }
      // }else if (sdeScheme == word("OrsteinTime") ) {
      //   dsigma_t = 0.5e0 * Foam::sqrt(dtau.value() * dtheta.value());
      //   forAll(U, i){
      //     dimensionedScalar dW = rndGen.scalarNormal();
      //     dimensionedScalar dWi = dW * dsigma_t;
      //     dr[i] = dW.value();
      //     if (i == 0){
      //         dF[i] = dWi.value();
      //     } else {
      //         dF[i] = rcorr * dF[i-1] + rroot * dWi.value();
      //     }
      //     velInit[i] = dF[i];
      //   }
      }else if (sdeScheme == word("OrsteinTimeN") ) {
        
        rcorr = Foam::exp(-dx*dtheta.value() / dnu);
        rroot = Foam::sqrt(1.0-rcorr*rcorr);

        dsigma_t = 0.5e0 * Foam::sqrt(dtau.value() * dtheta.value());
        dsigma_t /= Foam::sqrt(deps * dnu);
        forAll(U, i){
          // dimensionedScalar dW = rndGen.scalarNormal();
          dWs = rndGen.scalarNormal();
          dWi = dWs * dsigma_t;
          dr[i] = dWs;
          if (i == 0){
              dF[i] = dWi;
          } else {
              dF[i] = rcorr * dF[i-1] + rroot * dWi;
          }
          velInit[i] = dF[i];
        }
      } else if (sdeScheme == word("WhiteNoise")|| sdeScheme == word("WhiteNoiseTime")) {
        forAll(U, i){
          // dimensionedScalar dW = rndGen.scalarNormal();
          dWs = rndGen.scalarNormal();
          dWs /= Foam::sqrt(dx);
          dr[i] = dWs;
          velInit[i] = dWs;
        }
      } else if (sdeScheme == word("Poisson")) {

        Info << "Calculate RTN" << endl;
        // set the poisson distribution
        // std::default_random_engine generator;
        std::random_device rd;
        std::mt19937 gen(rd()); // clockTime.getTime()
        // std::mt19937 gen(clockTime.getTime()); // clockTime.getTime()
        std::poisson_distribution<int> distribution(1.0/dx);

        label n1 = label(1.0 * 10000/dx) * 10 ;
        label icnto = 0 ;
        label icnt = 0 ;
        scalar isgn = 1. ;
        // scalar dW = rndGen.scalarNormal();
        dWs = rndGen.scalarNormal();
        if (dWs < 0.5){
          isgn = -1.0;
         }
        Info << "First element "<< isgn << endl;
        for (label i = 0; i < n1; i++)
        {
          label np = distribution(gen);
          icnt += np;

          Info << "Gen " << np << " " << icnto << " " << icnt <<" " << isgn << endl;
          if (icnt > nel) {
            icnt = nel;
          }
          if ( np != 0) {
            for (label j=icnto;j<icnt; j++){
              dF[j] = isgn;
            }
          }
          icnto = icnt;
          if (isgn < 0){
            isgn = 1.0;
          }else {
            isgn = -1.0;
          }
          if (icnt >= nel) {
            break;
          }
        }
        
        forAll(U, i){
          velInit[i] = dF[i];
        }
      }
      // calc mean value of velInit and modify it to mean zero
      // scalar velm(0.0);
      // forAll(velInit, i){
      //     velm += velInit[i];
      // }
      // velm /= velInit.size();
      // forAll(velInit, i){
      //     velInit[i] -= velm;
      // }
    } else {

      Info << "Use existing stochastic field" << nl << endl;

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
    // std::ofstream file;
    // // file.open("res.dat",std::ofstream::out| std::ofstream::app);
    // file.open("res.dat",std::ofstream::out);
    // forAll(dF, i) {
    //   file << i<<","<< C[i][0] << ',' << velInit[i]<<","<<U[i][0] << ',' << dr[i]<<std::endl;
    // }
    // file.close();

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
        if (sdeScheme == word("OrsteinTime") ) {
          forAll(dF , i) {
            dFo[i] = dF[i];
          }

          dsigma_t = 0.5e0 * Foam::sqrt(dtau.value() * dtheta.value());
          //
          scalar rcorry(Foam::exp(-(runTime.deltaTValue())*dtau.value()));
          // timeo = runTime.deltaTValue();
          scalar rrooty(Foam::sqrt(1.0-rcorry*rcorry));
          forAll(U, i){
            // dimensionedScalar dW = rndGen.scalarNormal();
            dWs = rndGen.scalarNormal();
            // dWs = dnd(gen);
            dWi = dWs * dsigma_t;
            // dW = dW * dsigma * Foam::sqrt(dx);
            // dW /= Foam::sqrt(dx);
            // dr[i] = dW.value();
            if (i == 0){
    //           U[i][0] = dW.value() ;
                dF[i] = rcorry * dFo[i] + rrooty * dWi;
            } else {
                //U[i][0] = U[i-1][0] + dtheta.value()*(dmean.value() -U[i-1][0] ) * dx + dW.value();
                // U[i][0] = (1.0-dtheta.value()*dx) * U[i-1][0] + dW.value();
                dF[i] = rcorr * dF[i-1] + rcorry * dFo[i] - rcorr * rcorry * dFo[i-1] + rroot * rrooty * dWi;
            }
            velInit[i] = dF[i];

            U[i][0] = velInit[i] * barVel.value();
            U[i][0] -= vShift.value();

            // velInit[i] = dW.value();
          }
        } else if (sdeScheme == word("OrsteinTimeM") ) {
          forAll(dF , i) {
            dFo[i] = dF[i];
          }

          dsigma_t = 0.5e0 * Foam::sqrt(dtau.value() * dtheta.value());
          //
          scalar rcorry(Foam::exp(-(runTime.deltaTValue())*dtau.value()));
          // timeo = runTime.deltaTValue();
          scalar rrooty(Foam::sqrt(1.0-rcorry*rcorry));

          vel_mean = vel* 0.1 + vel_meano * 0.9;
          forAll(U, i){
            // dimensionedScalar dW = rndGen.scalarNormal();
            dWs = rndGen.scalarNormal();
            dWi = dWs * dsigma_t;
            // dW = dW * dsigma * Foam::sqrt(dx);
            // dW /= Foam::sqrt(dx);
            dr[i] = dWs;
            if (i == 0){
    //           U[i][0] = dW.value() ;
                dF[i] = rcorry * dFo[i] + rrooty * dWi;
            } else {
                //U[i][0] = U[i-1][0] + dtheta.value()*(dmean.value() -U[i-1][0] ) * dx + dW.value();
                // U[i][0] = (1.0-dtheta.value()*dx) * U[i-1][0] + dW.value();
                dF[i] = rcorr * dF[i-1] + rcorry * dFo[i] - rcorr * rcorry * dFo[i-1] + rroot * rrooty * dWi;
            }
            velInit[i] = dF[i];

            U[i][0] = velInit[i] * barVel.value();
            U[i][0] -= vShift.value();
            U[i][0] -= vel_mean.value() ;

            // velInit[i] = dW.value();
          }
          vel_meano = vel;
        } else if (sdeScheme == word("OrsteinTimeN") ) {
          forAll(dF , i) {
            dFo[i] = dF[i];
          }

          dsigma_t = 0.5e0 * Foam::sqrt(dtau.value() * dtheta.value());
          dsigma_t /= Foam::sqrt(deps * dnu);
          //
          scalar rcorry(Foam::exp(-(runTime.deltaTValue())*dtau.value() / deps));
          // timeo = runTime.deltaTValue();
          scalar rrooty(Foam::sqrt(1.0-rcorry*rcorry));
          forAll(U, i){
            // dimensionedScalar dW = rndGen.scalarNormal();
            dWs = rndGen.scalarNormal();
            dWi = dWs * dsigma_t;
            // dW = dW * dsigma * Foam::sqrt(dx);
            // dW /= Foam::sqrt(dx);
            dr[i] = dWs ;
            if (i == 0){
    //           U[i][0] = dW.value() ;
                dF[i] = rcorry * dFo[i] + rrooty * dWi;
            } else {
                //U[i][0] = U[i-1][0] + dtheta.value()*(dmean.value() -U[i-1][0] ) * dx + dW.value();
                // U[i][0] = (1.0-dtheta.value()*dx) * U[i-1][0] + dW.value();
                dF[i] = rcorr * dF[i-1] + rcorry * dFo[i] - rcorr * rcorry * dFo[i-1] + rroot * rrooty * dWi;
            }
            velInit[i] = dF[i];

            U[i][0] = velInit[i] * barVel.value();
            U[i][0] -= vShift.value();

            // velInit[i] = dW.value();
          }
        } else if (sdeScheme == word("WhiteNoiseTime")) {
          scalar dvar= Foam::sqrt(dx * runTime.deltaTValue());
          forAll(U, i){
            dimensionedScalar dW = rndGen.scalarNormal();
            dWs = rndGen.scalarNormal();
            dWs /= dvar;
            // dW /= Foam::sqrt(dx);
            // dW /= Foam::sqrt(runTime.deltaTValue());
            velInit[i] = dWs;
            U[i][0] = velInit[i] * barVel.value();
            U[i][0] -= vShift.value();
          }
        }

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
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              - fvm::div(phi, T)
              - fvm::laplacian(DT, T)
              + fvm::Sp(2*DK*T + fvc::div(U) - DK, T)
              ==
                // fvm::Sp(DK - DK*T - fvc::div(U), T)
                // fvm::Sp(DK - 2*DK*T - fvc::div(U), T)
            //   + fvc::Su(DK*T*T,T)
                fvc::Su(DK*T*T,T)
            );

            TEqn.relax();
            fvOptions.constrain(TEqn);
            // TEqn.solve();
            scalar TInitialResidual = TEqn.solve().finalResidual();
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

            // // cut-off the upper and lower bounds to 1 and 0
            // forAll(T,i) {
            //   T[i] = Foam::max(Foam::min(T[i] , 1) , 0);
            // }

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
        scalar flamePos(-1.0);
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
        // if (flamePos < 0.0) {
        //   flamePos  = C[i].component(0);
        // }
        meanVel /= mesh.C().size();
        Info << "ForCorr" ;
        for (int i=1;i<8;i++){
          // Info << " "<< U[i*100][0];
          printf(" %g",U[i*100][0]);
        }
        Info << endl;
        Info<< nl << "Velocity = " << vel.value()<< " " << velRes.value() << " " << tvel.value() << " " << flamePos << " " << Foam::gMax(T.internalField()) << " " << Foam::gMax(xi.internalField())<< " " << meanVel << nl << endl;      
        // Info<< "phi = " << phi.value()<< " " << velRes.value() << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
