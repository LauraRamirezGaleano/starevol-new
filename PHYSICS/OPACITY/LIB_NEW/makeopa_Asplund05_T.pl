#!/usr/bin/perl -w 

$DIRBASE = "/home/palacios/CALCULS/EVOL/PHYSICS/OPACITY/LIB_NEW/LTOPA";

$indchem = 0;

@files = ("ags04.0.1.tron",
	  "ags04.0.08.tron",
	  "ags04.0.06.tron",
	  "ags04.0.05.tron",
	  "ags04.0.04.tron",
	  "ags04.0.03.tron",
	  "ags04.0.02.tron",
	  "ags04.0.01.tron",
	  "ags04.0.004.tron",
	  "ags04.0.002.tron",
	  "ags04.0.001.tron",
	  "ags04.0.0003.tron",
	  "ags04.0.0001.tron",
	  "ags04.0.00003.tron",
	  "ags04.0.00001.tron",
	  "ags04.0.0.tron");
@files = (@files,"ags04.1.1.tron",
	  "ags04.1.08.tron",
	  "ags04.1.06.tron",
	  "ags04.1.05.tron",
	  "ags04.1.04.tron",
	  "ags04.1.03.tron",
	  "ags04.1.02.tron",
	  "ags04.1.01.tron",
	  "ags04.1.004.tron",
	  "ags04.1.002.tron",
	  "ags04.1.001.tron",
	  "ags04.1.0003.tron",
	  "ags04.1.0001.tron",
	  "ags04.1.00003.tron",
	  "ags04.1.00001.tron",
	  "ags04.1.0.tron");
@files = (@files,"ags04.2.1.tron",
	  "ags04.2.08.tron",
	  "ags04.2.06.tron",
	  "ags04.2.05.tron",
	  "ags04.2.04.tron",
	  "ags04.2.03.tron",
	  "ags04.2.02.tron",
	  "ags04.2.01.tron",
	  "ags04.2.004.tron",
	  "ags04.2.002.tron",
	  "ags04.2.001.tron",
	  "ags04.2.0003.tron",
	  "ags04.2.0001.tron",
	  "ags04.2.00003.tron",
	  "ags04.2.00001.tron",
	  "ags04.2.0.tron");
@files = (@files,"ags04.35.1.tron",
	  "ags04.35.08.tron",
	  "ags04.35.06.tron",
	  "ags04.35.05.tron",
	  "ags04.35.04.tron",
	  "ags04.35.03.tron",
	  "ags04.35.02.tron",
	  "ags04.35.01.tron",
	  "ags04.35.004.tron",
	  "ags04.35.002.tron",
	  "ags04.35.001.tron",
	  "ags04.35.0003.tron",
	  "ags04.35.0001.tron",
	  "ags04.35.00003.tron",
	  "ags04.35.00001.tron",
	  "ags04.35.0.tron");
@files = (@files,"ags04.5.1.tron",
	  "ags04.5.08.tron",
	  "ags04.5.06.tron",
	  "ags04.5.05.tron",
	  "ags04.5.04.tron",
	  "ags04.5.03.tron",
	  "ags04.5.02.tron",
	  "ags04.5.01.tron",
	  "ags04.5.004.tron",
	  "ags04.5.002.tron",
	  "ags04.5.001.tron",
	  "ags04.5.0003.tron",
	  "ags04.5.0001.tron",
	  "ags04.5.00003.tron",
	  "ags04.5.00001.tron",
	  "ags04.5.0.tron");
@files = (@files,"ags04.7.1.tron",
	  "ags04.7.08.tron",
	  "ags04.7.06.tron",
	  "ags04.7.05.tron",
	  "ags04.7.04.tron",
	  "ags04.7.03.tron",
	  "ags04.7.02.tron",
	  "ags04.7.01.tron",
	  "ags04.7.004.tron",
	  "ags04.7.002.tron",
	  "ags04.7.001.tron",
	  "ags04.7.0003.tron",
	  "ags04.7.0001.tron",
	  "ags04.7.00003.tron",
	  "ags04.7.00001.tron",
	  "ags04.7.0.tron");
@files = (@files,"ags04.8.1.tron",
	  "ags04.8.08.tron",
	  "ags04.8.06.tron",
	  "ags04.8.05.tron",
	  "ags04.8.04.tron",
	  "ags04.8.03.tron",
	  "ags04.8.02.tron",
	  "ags04.8.01.tron",
	  "ags04.8.004.tron",
	  "ags04.8.002.tron",
	  "ags04.8.001.tron",
	  "ags04.8.0003.tron",
	  "ags04.8.0001.tron",
	  "ags04.8.00003.tron",
	  "ags04.8.00001.tron",
	  "ags04.8.0.tron");
@files = (@files,"ags04.9.1.tron",
	  "ags04.9.08.tron",
	  "ags04.9.06.tron",
	  "ags04.9.05.tron",
	  "ags04.9.04.tron",
	  "ags04.9.03.tron",
	  "ags04.9.02.tron",
	  "ags04.9.01.tron",
	  "ags04.9.004.tron",
	  "ags04.9.002.tron",
	  "ags04.9.001.tron",
	  "ags04.9.0003.tron",
	  "ags04.9.0001.tron",
	  "ags04.9.00003.tron",
	  "ags04.9.00001.tron",
	  "ags04.9.0.tron");
@files = (@files,"ags04.92.08.tron");
@files = (@files,"ags04.94.06.tron");
@files = (@files,"ags04.95.0.tron",
	  "ags04.95.05.tron",
	  "ags04.95.04.tron",
	  "ags04.95.03.tron",
	  "ags04.95.02.tron",
	  "ags04.95.01.tron",
	  "ags04.95.004.tron",
	  "ags04.95.002.tron",
	  "ags04.95.001.tron",
	  "ags04.95.0003.tron",
	  "ags04.95.0001.tron",
	  "ags04.95.00003.tron",
	  "ags04.95.00001.tron");
@files = (@files,"ags04.96.04.tron");
@files = (@files,"ags04.97.03.tron");
@files = (@files,"ags04.98.02.tron");
@files = (@files,"ags04.99.01.tron");
@files = (@files,"ags04.996.004.tron");
@files = (@files,"ags04.998.002.tron");
@files = (@files,"ags04.999.001.tron");
@files = (@files,"ags04.9997.0003.tron");
@files = (@files,"ags04.9999.0001.tron");
@files = (@files,"ags04.99997.00003.tron");
@files = (@files,"ags04.99999.00001.tron");
@files = (@files,"ags04.one.tron");


print ("      BLOCK DATA LTHOPA

************************************************************************
*                     Tables of radiative opacities                    *
*                   at low temperature (T < 35000 K)                   *
*                          Ferguson (2005)                             *
************************************************************************

      implicit none

      double precision rklth,tklth,opaclth
      integer km,kt

      common /lthopac/ rklth(19),tklth(85),opaclth(19,85,155)

      data (rklth(km),km = 1,19)/ -8.000d0,-7.500d0,-7.000d0,-6.500d0,
     & -6.000d0,-5.500d0,-5.000d0,-4.500d0,-4.000d0,-3.500d0,-3.000d0,
     & -2.500d0,-2.000d0,-1.500d0,-1.000d0,-0.500d0, 0.000d0, 0.500d0,
     &  1.000d0 /
      data (tklth(kt),kt = 1,85)/3.162278d4,2.818383d4,2.511886d4,
     &2.238721d4,1.995262d4,1.778279d4,1.584893d4,1.412538d4,1.258925d4,
     &1.122018d4,1.000000d4,8.912509d3,7.943282d3,7.079458d3,6.309573d3,
     &5.623413d3,5.011872d3,4.466836d3,3.981072d3,3.548134d3,3.162278d3,
     &3.090295d3,3.019952d3,2.951209d3,2.884032d3,2.818383d3,2.754229d3,
     &2.691535d3,2.630268d3,2.570396d3,2.511886d3,2.454709d3,2.398833d3,
     &2.344229d3,2.290868d3,2.238721d3,2.187762d3,2.137962d3,2.089296d3,
     &2.041738d3,1.995262d3,1.949845d3,1.905461d3,1.862087d3,1.819701d3,
     &1.778279d3,1.737801d3,1.698244d3,1.659587d3,1.621810d3,1.584893d3,
     &1.548817d3,1.513561d3,1.479108d3,1.445440d3,1.412538d3,1.380384d3,
     &1.348963d3,1.318257d3,1.288250d3,1.258925d3,1.230269d3,1.202264d3,
     &1.174898d3,1.148154d3,1.122018d3,1.096478d3,1.071519d3,1.047129d3,
     &1.023293d3,1.000000d3,9.772372d2,9.549926d2,9.332543d2,9.120108d2,
     &8.912509d2,8.709636d2,8.511380d2,8.317638d2,8.128305d2,7.943282d2,
     &7.079458d2,6.309573d2,5.623413d2,5.011872d2 \/ \n\n");

foreach (@files)
{
    chomp;
    open (OPA,"$DIRBASE/$_");
    $indchem++;
    $indT = 0;
    foreach (<OPA>)
    {
	if (/^G/)
	{
	    print ("*\n*\t Table $_*\n");
	}
	elsif (/^[234]\./)
	{
	    s/-1/ -1/g;
	    @val = split(/\s+/);
	    $indT++;
	    print ("*\t\t\t log(T) = $val[0]\n");
	    print ("      data (opaclth(km,$indT,$indchem),km = 1,19)\/");
	    print ("$val[1]d0, $val[2]d0,\n");
	    print ("     &\t$val[3]d0,$val[4]d0, $val[5]d0, $val[6]d0, $val[7]d0, $val[8]d0,\n");
	    print ("     &\t$val[9]d0, $val[10]d0, $val[11]d0, $val[12]d0, $val[13]d0, $val[14]d0,\n");
	    print ("     &\t$val[15]d0, $val[16]d0, $val[17]d0, $val[18]d0, $val[19]d0\/\n");
	}	
    }
    close(OPA);    
}
print ("      end\n");


