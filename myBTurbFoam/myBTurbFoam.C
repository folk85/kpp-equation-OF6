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
#include <complex>
#include <iostream>
#include <valarray>
#include "complexFields.H"

const double PI = 3.141592653589793238460;
 
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

#define NMODES 3

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;
 
    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];
 
    // conquer
    fft(even);
    fft(odd);
 
    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

// inverse fft (in-place)
void ifft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);
 
    // forward fft
    fft( x );
 
    // conjugate the complex numbers again
    x = x.apply(std::conj);
 
    // scale the numbers
    x /= x.size();
}

scalar get_ee(scalar kk, scalar theta){
  return theta / 2.0 / 3.1415 / (kk * kk + theta * theta);
}


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
    Info << "nmodes_read = " << nmodes_read << endl ;
    Info << "mean_period_iters = " << mean_period_iters << endl ;
    Info << "ampl_forcing = " << ampl_forcing << endl ;

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
    scalar dnel = scalar(nel);
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
      if ((sdeScheme == word("Orstein"))||(sdeScheme == word("OrsteinTime"))||(sdeScheme == word("OrsteinTimeF1") ) ) {
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

    // Set the field 
    CArray ufft(C.size());
    const double imag = 0.0e0;
    forAll(U, i) {
      std::complex<double> ctmp;
      ctmp.real(U[i][0]);
      ctmp.imag(imag);
      ufft[i] = ctmp;
      // ufft[i].imag = 0.0;
    }
    fft(ufft);

    List<scalar> alpha(3);
    alpha[0] = 1.0e-1;
    alpha[1] = 1.0e-3;
    alpha[2] = 100;
    
    complexField ek_c0(C.size());
    complexField ek_c1(C.size());
    complexField ek_c2(C.size());
    complexField ek_c3(C.size());
    CArray r_c1(C.size());
    CArray r_c2(C.size());
    CArray r_c3(C.size());
    forAll(ek_c0, i) {
      // ek_c[i] = ufft[i] * std::conj(ufft[i]);
      ek_c0[i] = ufft[i].real() * ufft[i].real() + ufft[i].imag() * ufft[i].imag();
      ek_c0[i] /= dnel;
      ek_c1[i] = ek_c0[i];
      ek_c2[i] = ek_c0[i];
      ek_c3[i] = ek_c0[i];
      // r_c1[i].real(ek_c1[i].Re());
      // r_c2[i].real(ek_c2[i].Re());
      // r_c3[i].real(ek_c3[i].Re());
      // r_c1[i].imag(imag);
      // r_c2[i].imag(imag);
      // r_c3[i].imag(imag);
    }
    // ifft(r_c1);
    // ifft(r_c2);
    // ifft(r_c3);
    // forAll(r_c1, i) {
    //   r_c1[i] /= r_c1[0];
    //   r_c2[i] /= r_c2[0];
    //   r_c3[i] /= r_c3[0];
    // }

    #include "setInitialDeltaT.H"

    // set forcing coefficients for Burgers equation

    Info << "Set forcing sources" << endl;

    std::random_device rdt;
    std::mt19937 gen_ts(rdt()); // clockTime.getTime()

    // std::poisson_distribution<int> poisson_delta(10.0 / dx);
    std::poisson_distribution<int> poisson_delta(mean_period_iters);

    scalar dfreq = 2.0 * 3.1415 / dx / C.size();
    List<scalar> kfreq(nmodes_read);
    for (int i = 0; i < nmodes_read; i++)
    {
      kfreq[i] = dfreq * (i + 1);
    }
    
    List<scalar> qamp(nmodes_read);
    scalar qamptmp = dtheta.value() * barVel.value();
    for (int i = 0; i < nmodes_read; i++)
    {
      qamp[i] = qamptmp * 0.5 / 3.1415 / (kfreq[i] * kfreq[i] + qamptmp * qamptmp);
      qamp[i] = Foam::sqrt(qamp[i] * dfreq);
    }
    List<scalar> modesa(nmodes_read);
    List<scalar> modesb(nmodes_read);
    scalar time_change = runTime.value();
    scalar dtime_delta = dx * dx * 0.8 ;
    label nsteps = poisson_delta(gen_ts);
    scalar maxa(-1.0);
    scalar maxb(-1.0);
    
    scalar magnitude = ampl_forcing;
    
    time_change += dtime_delta * nsteps;
    printf("Initialize forcing modes. Get N_iters: %d , Next update after: %g\n", nsteps, time_change);

    for (int i = 0; i < nmodes_read; i++)
    {
      modesa[i] = rndGen.scalarNormal();
      modesb[i] = rndGen.scalarNormal();
      // if (abs(modesa[i]) > maxa)
      // maxa = abs(modesa[i]);
      // if (abs(modesb[i]) > maxb)
      // maxb = abs(modesb[i]);
    }
    // for (int i = 0; i < nmodes_read; i++)
    // {
    //   modesa[i] /= maxa;
    //   modesb[i] /= maxb;
    // }
    Info << "volVectorField forcing" << endl;
    for (int i = 0; i < C.size(); i++)
    {
      scalar dtmp = 0.0e0;
      for (int j = 0; j < nmodes_read; j++)
      {
        dtmp += qamp[j] * (modesa[j] * Foam::cos(kfreq[j] * C[i].component(0)) + modesb[j] * Foam::sin(kfreq[j] * C[i].component(0)));
      }
      // printf("forcing[%d] = %g  cos = %g   sin = %g\n", i, dtmp, Foam::cos(kfreq[i] * C[i].component(0)), Foam::sin(kfreq[i] * C[i].component(0)));
      forcing[i][0] = dtmp * magnitude;
      // forcing[i][0] = 0.0;
      forcing[i][1] = 0.0;
      forcing[i][2] = 0.0;
    }
    
    
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

        // Update random variables for forcing
        if (runTime.value() >= time_change){

          nsteps = poisson_delta(gen_ts);
          // Info << "update forcing modes " << nsteps << endl;
          time_change = runTime.value() + dtime_delta * nsteps;
          printf("Update forcing modes. Get N_iters: %d , Next update after: %g\n", nsteps, time_change);
          maxa = -1.0;
          maxb = -1.0;
          for (int i = 0; i < nmodes_read; i++)
          {
            modesa[i] = rndGen.scalarNormal();
            modesb[i] = rndGen.scalarNormal();
            // if (abs(modesa[i]) > maxa)
            //   maxa = abs(modesa[i]);
            // if (abs(modesb[i]) > maxb)
            //   maxb = abs(modesb[i]);
          }
          // for (int i = 0; i < nmodes_read; i++)
          // {
          //   modesa[i] /= maxa;
          //   modesb[i] /= maxb;
          // }
          for (int i = 0; i < C.size(); i++)
          {
            scalar dtmp = 0.0e0;
            for (int j = 0; j < nmodes_read; j++)
            {
              dtmp += qamp[j] * (modesa[j] * Foam::cos(kfreq[j] * C[i].component(0)) + modesb[j] * Foam::sin(kfreq[j] * C[i].component(0)));
            }
            forcing[i][0] = dtmp * magnitude;
            // forcing[i][0] = 0.0;
            forcing[i][1] = 0.0;
            forcing[i][2] = 0.0;
          }
          Info << "End updating forcing modes" << endl;
        }

//             // Update velocity
        if ((sdeScheme == word("OrsteinTime") )|| (sdeScheme == word("OrsteinTimeF1") )) {
          forAll(dF , i) {
            dFo[i] = dF[i];
          }

          dsigma_t = 0.5e0 * Foam::sqrt(dtau.value() * dtheta.value());
          scalar rcorry(Foam::exp(-(runTime.deltaTValue())*dtau.value()));
          if ( sdeScheme == word("OrsteinTimeF1") ) {
            scalar dtau_new = dtheta.value() * 0.25e0 * barVel.value() * dtheta.value() * barVel.value() * dtheta.value();
            dsigma_t = 0.5e0 * Foam::sqrt(dtau_new * dtheta.value());
            rcorry = Foam::exp(-(runTime.deltaTValue()) * dtau_new);
          }         
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
      Vector<scalar> UInitialResidual;
      for (int i = 0; i < 5; i++) {
        fvVectorMatrix UEqn
        (
       fvm::ddt(U)
           + 0.5 * fvm::div(phi, U)
           - fvm::laplacian(nu, U)
           ==
          fvc::Su(forcing, U)
        );

        UEqn.relax();
    // fvOptions.constrain(UEqn);
        UInitialResidual = UEqn.solve().finalResidual();
        // phi = fvc::flux(U);
        phi = fvc::interpolate(U) & mesh.Sf();
      }
        // fvOptions.correct(UEqn);


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

            fvScalarMatrix xiEqn
            (
                fvm::ddt(xi)
              - fvm::div(phi, xi)
              - fvm::laplacian(DT, xi)
              + fvm::Sp(2*DK*xi + fvc::div(U) - DK, xi)
              ==
                // fvm::Sp(DK - DK*T - fvc::div(U), T)
                // fvm::Sp(DK - 2*DK*T - fvc::div(U), T)
            //   + fvc::Su(DK*T*T,T)
                fvc::Su(DK*xi*xi,xi)
            );

            xiEqn.relax();
            fvOptions.constrain(xiEqn);
            scalar xiInitialResidual = xiEqn.solve().finalResidual();
            fvOptions.correct(xi);

            Info<< "Residuals = " << TInitialResidual << " " << xiInitialResidual << " Cnt = "<< icount  << endl;

            if (TInitialResidual<1.0e-5 && xiInitialResidual < 1.0e-5 && UInitialResidual[0] < 1.0e-5) {
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
        scalar meanVel2(0.0);
        scalar stdVel(0.0);
        forAll(mesh.C(),i) {
            meanVel += U[i][0];
            meanVel2 += magSqr(U[i]);
        }
        forAll(mesh.C(),i) {
            if (T[i] >= 0.50e0) {
                flamePos = C[i].component(0);
                break;
            }
        }
        flamePos = gSum(T) * dx;
        scalar flamePos2 = gSum(xi) * dx;
        // if (flamePos < 0.0) {
        //   flamePos  = C[i].component(0);
        // }
        meanVel /= mesh.C().size();
        meanVel2 /= mesh.C().size();
        forAll(mesh.C(),i) {
            stdVel += (U[i][0] - meanVel) * (U[i][0] - meanVel);
        }
        stdVel /= mesh.C().size();
        stdVel = Foam::sqrt(stdVel);
        Info << "ForCorr" ;
        for (int i=1;i<8;i++){
          // Info << " "<< U[i*100][0];
          printf(" %g",U[i*100][0]);
        }
        Info << endl;
        //
        // Calc Spectral energy of 
        forAll(U, i) {
          std::complex<double> ctmp;
          ctmp.real(U[i][0]);
          ctmp.imag(imag);
          ufft[i] = ctmp;
          // ufft[i].imag = 0.0;
        }
        fft(ufft);
        scalar alpha_t = runTime.deltaTValue() / alpha[2];
        forAll(ek_c0, i) {
          // ek_c[i] = ufft[i] * std::conj(ufft[i]);
          ek_c0[i] = ufft[i].real() * ufft[i].real() + ufft[i].imag() * ufft[i].imag();
          ek_c0[i] /= dnel;
          ek_c1[i] = ek_c1[i] * (1.0 - alpha[0]) + alpha[0] * ek_c0[i];
          ek_c2[i] = ek_c2[i] * (1.0 - alpha[1]) + alpha[1] * ek_c0[i];
          ek_c3[i] = ek_c3[i] * (1.0 - alpha_t) + alpha_t * ek_c0[i];
          Ek[i][0] = ek_c1[i].Re();
          Ek[i][1] = ek_c2[i].Re();
          Ek[i][2] = ek_c3[i].Re();
          r_c1[i].real(ek_c1[i].Re());
          r_c2[i].real(ek_c2[i].Re());
          r_c3[i].real(ek_c3[i].Re());
          r_c1[i].imag(imag);
          r_c2[i].imag(imag);
          r_c3[i].imag(imag);
        }
        ifft(r_c1);
        ifft(r_c2);
        ifft(r_c3);
        forAll(r_c1, i) {
          r_c1[i] /= r_c1[0];
          r_c2[i] /= r_c2[0];
          r_c3[i] /= r_c3[0];
          Rx[i][0] = (r_c1[i] / r_c1[0]).real();
          Rx[i][1] = (r_c2[i] / r_c2[0]).real();
          Rx[i][2] = (r_c3[i] / r_c3[0]).real();
        }
        //
        // Calc L = \int_{0}^{Lx/2} Rx dx
        //
        std::vector<double> lx(6);
        forAll(lx, i) {
          lx[i] = 0.0e0;
        }
        for (label i = 0; i < nel / 2; i++)
        {
          lx[0] += dx * Rx[i][0];
          lx[1] += dx * Rx[i][1];
          lx[2] += dx * Rx[i][2];
          if (i < nel / 4) {
            lx[3] += dx * Rx[i][0];
            lx[4] += dx * Rx[i][1];
            lx[5] += dx * Rx[i][2];
          }
        }
        
        printf("Correlation L ");
        forAll(lx, i) {
          printf("%g ", lx[i]);
        }
        printf("%g ", Ek[0][0]);
        printf("%g ", Ek[0][1]);
        printf("%g ", Ek[0][2]);
        printf("%g ", r_c1[0].real());
        printf("%g ", r_c2[0].real());
        printf("%g\n", r_c3[0].real());

        printf("Get the velocity parameters Mean %g   MeanSqr %g and stdVel %g \n", meanVel, meanVel2, stdVel);
        // Info<< nl << "Velocity = " << vel.value()<< " " << velRes.value() << " " << tvel.value() << " " << flamePos << " " << Foam::gMax(T.internalField()) << " " << Foam::gMax(xi.internalField())<< " " << meanVel << " " << flamePos2 << " " << meanVel2 << nl << endl;      
        // Info<< "phi = " << phi.value()<< " " << velRes.value() << nl << endl;

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
