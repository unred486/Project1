C**********************************************************************C
C                                                                      C
C   soot.h (03/23/99)                                                  C
C   Include file for soot code                                         C
C                                                                      C
C**********************************************************************C

	integer HI_MOMENT, MAX_MOMENTS
	parameter(HI_MOMENT = 5, MAX_MOMENTS = HI_MOMENT + 1)

C	Highest and lowest fractional moments to be calculated
	integer LO_FRAC, HI_FRAC
	parameter(LO_FRAC = -4, HI_FRAC = 6 * MAX_MOMENTS + 1)

	double precision PI
	parameter(PI = 3.141592654)
!	parameter(PI = 3.14159)

	double precision ONE_SIXTH
	parameter(ONE_SIXTH = 1.0D0 / 6.0D0)

	double precision ONE_THIRD
	parameter(ONE_THIRD = 1.0D0 / 3.0D0)

	double precision TWO_THIRDS
	parameter(TWO_THIRDS = 2.0D0 / 3.0D0)

	double precision AVOGADRO
	parameter(AVOGADRO = 6.022136D+23)

	double precision BOLTZMANN
	parameter(BOLTZMANN = 1.3807D-16)

	double precision GAS_CONST
	parameter(GAS_CONST = 8.314D-03)

C	Atomic mass unit (g)
	double precision AMU
	parameter(AMU = 1.67D-24)

C	Mass of carbon atom (g)
	double precision C_MASS
	parameter(C_MASS = 12.0D0 * AMU)

C	Size of the benzene ring in (cm)
C	2.42 = 1.395*sqrt(3), 1.395 A is the C-C in aromatics
	double precision D1_PAH
	parameter(D1_PAH = 2.42D-8)

C	Density of soot in g/cm^3
	double precision RHO_SOOT
	parameter(RHO_SOOT = 1.8)

C	Factor used in coagulation (contimuum regime - viscosity term)
	double precision CK
	parameter (CK = 2.0 * BOLTZMANN / 3.0)

C	Nominal # of sites per cm^2
	double precision CHI
	parameter(CHI = 2.3D+15)

C	Mass of OH (g)
	double precision OH_MASS
	parameter(OH_MASS = 17.0 * AMU)

C	Lower limit of M0 during oxidation (arbitrary)
	double precision M0_OXID_LIMIT
	parameter(M0_OXID_LIMIT = 1.0D0)

C	Constant used to calculate M0 rate during oxidation (empirical)
	double precision M0_OXID_CONST
	parameter(M0_OXID_CONST = 5.0D-11)

C	Common block for inputs to initSoot

	integer nMoments, jPAH, jC2H2, jCO, jH, jH2, jH2O,
     &	                                    jO2, jOH, jM0
	logical fmOnly
	double precision cutoff
	double precision CPerPAH
	double precision diamPAH
	logical first
	common /INIT/ nMoments, jPAH, jC2H2, jCO, jH, jH2, jH2O,
     &	              jO2, jOH, jM0, fmOnly, first, cutoff, diamPAH,
     &	              CPerPAH

C	Common block for coefficients, etc. calculated in initSoot

	double precision CD1, cbNucl, CB0, CBCOND, CBOH, SURFC, CA,
     &	                 powPAH(3,5), x(MAX_MOMENTS),
     &	                 prime(MAX_MOMENTS,MAX_MOMENTS)
	double precision CPerDimer, binom(0:MAX_MOMENTS,0:MAX_MOMENTS)
	common / COEFFS / CD1, cbNucl, CB0, CBCOND, CBOH, SURFC, CA,
     &	                  powPAH, x, prime, CPerDimer, binom
