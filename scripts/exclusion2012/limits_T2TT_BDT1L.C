float getUpperLimit_T2TT_BDT1L( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 146.7;
     if(seff >= 0.05 && seff < 0.10) ul = 147.1;
     if(seff >= 0.10 && seff < 0.15) ul = 150.1;
     if(seff >= 0.15 && seff < 0.20) ul = 152.6;
     if(seff >= 0.20 && seff < 0.25) ul = 156.4;
     if(seff >= 0.25 && seff < 0.30) ul = 159.9;
     if(seff >= 0.30 && seff < 0.35) ul = 166.5;
     if(seff >= 0.35 && seff < 0.40) ul = 171.6;
     if(seff >= 0.40 && seff < 0.45) ul = 178.3;
     if(seff >= 0.45 && seff < 0.50) ul = 184.7;
     if(seff >= 0.50 && seff < 0.55) ul = 192.7;
     return ul;
}
float getExpectedUpperLimit_T2TT_BDT1L( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 164.5;
     if(seff >= 0.05 && seff < 0.10) ul = 163.1;
     if(seff >= 0.10 && seff < 0.15) ul = 165.9;
     if(seff >= 0.15 && seff < 0.20) ul = 169.2;
     if(seff >= 0.20 && seff < 0.25) ul = 175.0;
     if(seff >= 0.25 && seff < 0.30) ul = 180.3;
     if(seff >= 0.30 && seff < 0.35) ul = 187.5;
     if(seff >= 0.35 && seff < 0.40) ul = 193.7;
     if(seff >= 0.40 && seff < 0.45) ul = 201.9;
     if(seff >= 0.45 && seff < 0.50) ul = 208.5;
     if(seff >= 0.50 && seff < 0.55) ul = 216.4;
     return ul;
}
float getExpectedP1UpperLimit_T2TT_BDT1L( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 219.6;
     if(seff >= 0.05 && seff < 0.10) ul = 222.0;
     if(seff >= 0.10 && seff < 0.15) ul = 228.8;
     if(seff >= 0.15 && seff < 0.20) ul = 237.6;
     if(seff >= 0.20 && seff < 0.25) ul = 249.3;
     if(seff >= 0.25 && seff < 0.30) ul = 263.4;
     if(seff >= 0.30 && seff < 0.35) ul = 279.9;
     if(seff >= 0.35 && seff < 0.40) ul = 294.8;
     if(seff >= 0.40 && seff < 0.45) ul = 312.7;
     if(seff >= 0.45 && seff < 0.50) ul = 333.8;
     if(seff >= 0.50 && seff < 0.55) ul = 356.2;
     return ul;
}
float getExpectedM1UpperLimit_T2TT_BDT1L( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 120.5;
     if(seff >= 0.05 && seff < 0.10) ul = 120.4;
     if(seff >= 0.10 && seff < 0.15) ul = 121.9;
     if(seff >= 0.15 && seff < 0.20) ul = 122.5;
     if(seff >= 0.20 && seff < 0.25) ul = 124.0;
     if(seff >= 0.25 && seff < 0.30) ul = 126.0;
     if(seff >= 0.30 && seff < 0.35) ul = 129.6;
     if(seff >= 0.35 && seff < 0.40) ul = 132.7;
     if(seff >= 0.40 && seff < 0.45) ul = 136.2;
     if(seff >= 0.45 && seff < 0.50) ul = 139.4;
     if(seff >= 0.50 && seff < 0.55) ul = 142.9;
     return ul;
}
