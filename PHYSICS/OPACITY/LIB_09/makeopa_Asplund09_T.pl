#!/usr/bin/perl -w 

$DIRBASE = "/home/ana/CALCULS/EVOL/PHYSICS/OPACITY/LIB_09/LTOPA";

$indchem = 0;

@files = ("A09photo.0.1.tron",
	  "A09photo.0.08.tron",
	  "A09photo.0.06.tron",
	  "A09photo.0.05.tron",
	  "A09photo.0.04.tron",
	  "A09photo.0.03.tron",
	  "A09photo.0.02.tron",
	  "A09photo.0.01.tron",
	  "A09photo.0.004.tron",
	  "A09photo.0.002.tron",
	  "A09photo.0.001.tron",
	  "A09photo.0.0003.tron",
	  "A09photo.0.0001.tron",
	  "A09photo.0.00003.tron",
	  "A09photo.0.00001.tron",
	  "A09photo.0.0.tron");
@files = (@files,"A09photo.1.1.tron",
	  "A09photo.1.08.tron",
	  "A09photo.1.06.tron",
	  "A09photo.1.05.tron",
	  "A09photo.1.04.tron",
	  "A09photo.1.03.tron",
	  "A09photo.1.02.tron",
	  "A09photo.1.01.tron",
	  "A09photo.1.004.tron",
	  "A09photo.1.002.tron",
	  "A09photo.1.001.tron",
	  "A09photo.1.0003.tron",
	  "A09photo.1.0001.tron",
	  "A09photo.1.00003.tron",
	  "A09photo.1.00001.tron",
	  "A09photo.1.0.tron");
@files = (@files,"A09photo.2.1.tron",
	  "A09photo.2.08.tron",
	  "A09photo.2.06.tron",
	  "A09photo.2.05.tron",
	  "A09photo.2.04.tron",
	  "A09photo.2.03.tron",
	  "A09photo.2.02.tron",
	  "A09photo.2.01.tron",
	  "A09photo.2.004.tron",
	  "A09photo.2.002.tron",
	  "A09photo.2.001.tron",
	  "A09photo.2.0003.tron",
	  "A09photo.2.0001.tron",
	  "A09photo.2.00003.tron",
	  "A09photo.2.00001.tron",
	  "A09photo.2.0.tron");
@files = (@files,"A09photo.35.1.tron",
	  "A09photo.35.08.tron",
	  "A09photo.35.06.tron",
	  "A09photo.35.05.tron",
	  "A09photo.35.04.tron",
	  "A09photo.35.03.tron",
	  "A09photo.35.02.tron",
	  "A09photo.35.01.tron",
	  "A09photo.35.004.tron",
	  "A09photo.35.002.tron",
	  "A09photo.35.001.tron",
	  "A09photo.35.0003.tron",
	  "A09photo.35.0001.tron",
	  "A09photo.35.00003.tron",
	  "A09photo.35.00001.tron",
	  "A09photo.35.0.tron");
@files = (@files,"A09photo.5.1.tron",
	  "A09photo.5.08.tron",
	  "A09photo.5.06.tron",
	  "A09photo.5.05.tron",
	  "A09photo.5.04.tron",
	  "A09photo.5.03.tron",
	  "A09photo.5.02.tron",
	  "A09photo.5.01.tron",
	  "A09photo.5.004.tron",
	  "A09photo.5.002.tron",
	  "A09photo.5.001.tron",
	  "A09photo.5.0003.tron",
	  "A09photo.5.0001.tron",
	  "A09photo.5.00003.tron",
	  "A09photo.5.00001.tron",
	  "A09photo.5.0.tron");
@files = (@files,"A09photo.7.1.tron",
	  "A09photo.7.08.tron",
	  "A09photo.7.06.tron",
	  "A09photo.7.05.tron",
	  "A09photo.7.04.tron",
	  "A09photo.7.03.tron",
	  "A09photo.7.02.tron",
	  "A09photo.7.01.tron",
	  "A09photo.7.004.tron",
	  "A09photo.7.002.tron",
	  "A09photo.7.001.tron",
	  "A09photo.7.0003.tron",
	  "A09photo.7.0001.tron",
	  "A09photo.7.00003.tron",
	  "A09photo.7.00001.tron",
	  "A09photo.7.0.tron");
@files = (@files,"A09photo.8.1.tron",
	  "A09photo.8.08.tron",
	  "A09photo.8.06.tron",
	  "A09photo.8.05.tron",
	  "A09photo.8.04.tron",
	  "A09photo.8.03.tron",
	  "A09photo.8.02.tron",
	  "A09photo.8.01.tron",
	  "A09photo.8.004.tron",
	  "A09photo.8.002.tron",
	  "A09photo.8.001.tron",
	  "A09photo.8.0003.tron",
	  "A09photo.8.0001.tron",
	  "A09photo.8.00003.tron",
	  "A09photo.8.00001.tron",
	  "A09photo.8.0.tron");
@files = (@files,"A09photo.9.1.tron",
	  "A09photo.9.08.tron",
	  "A09photo.9.06.tron",
	  "A09photo.9.05.tron",
	  "A09photo.9.04.tron",
	  "A09photo.9.03.tron",
	  "A09photo.9.02.tron",
	  "A09photo.9.01.tron",
	  "A09photo.9.004.tron",
	  "A09photo.9.002.tron",
	  "A09photo.9.001.tron",
	  "A09photo.9.0003.tron",
	  "A09photo.9.0001.tron",
	  "A09photo.9.00003.tron",
	  "A09photo.9.00001.tron",
	  "A09photo.9.0.tron");
@files = (@files,"A09photo.92.08.tron");
@files = (@files,"A09photo.94.06.tron");
@files = (@files,"A09photo.95.0.tron",
	  "A09photo.95.05.tron",
	  "A09photo.95.04.tron",
	  "A09photo.95.03.tron",
	  "A09photo.95.02.tron",
	  "A09photo.95.01.tron",
	  "A09photo.95.004.tron",
	  "A09photo.95.002.tron",
	  "A09photo.95.001.tron",
	  "A09photo.95.0003.tron",
	  "A09photo.95.0001.tron",
	  "A09photo.95.00003.tron",
	  "A09photo.95.00001.tron");
@files = (@files,"A09photo.96.04.tron");
@files = (@files,"A09photo.97.03.tron");
@files = (@files,"A09photo.98.02.tron");
@files = (@files,"A09photo.99.01.tron");
@files = (@files,"A09photo.996.004.tron");
@files = (@files,"A09photo.998.002.tron");
@files = (@files,"A09photo.999.001.tron");
@files = (@files,"A09photo.9997.0003.tron");
@files = (@files,"A09photo.9999.0001.tron");
@files = (@files,"A09photo.99997.00003.tron");
@files = (@files,"A09photo.99999.00001.tron");
@files = (@files,"A09photo.10.0.tron");

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
	if (/^A/)
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


