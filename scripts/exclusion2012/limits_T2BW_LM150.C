float getUpperLimit_T2BW_LM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 95.6;
     if(seff >= 0.05 && seff < 0.10) ul = 96.4;
     if(seff >= 0.10 && seff < 0.15) ul = 97.9;
     if(seff >= 0.15 && seff < 0.20) ul = 99.5;
     if(seff >= 0.20 && seff < 0.25) ul = 101.4;
     if(seff >= 0.25 && seff < 0.30) ul = 103.6;
     if(seff >= 0.30 && seff < 0.35) ul = 106.6;
     if(seff >= 0.35 && seff < 0.40) ul = 109.5;
     if(seff >= 0.40 && seff < 0.45) ul = 112.7;
     if(seff >= 0.45 && seff < 0.50) ul = 116.3;
     if(seff >= 0.50 && seff < 0.55) ul = 120.3;
     return ul;
}
float getExpectedUpperLimit_T2BW_LM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 119.4;
     if(seff >= 0.05 && seff < 0.10) ul = 119.4;
     if(seff >= 0.10 && seff < 0.15) ul = 121.3;
     if(seff >= 0.15 && seff < 0.20) ul = 123.9;
     if(seff >= 0.20 && seff < 0.25) ul = 127.1;
     if(seff >= 0.25 && seff < 0.30) ul = 131.3;
     if(seff >= 0.30 && seff < 0.35) ul = 135.8;
     if(seff >= 0.35 && seff < 0.40) ul = 140.2;
     if(seff >= 0.40 && seff < 0.45) ul = 145.1;
     if(seff >= 0.45 && seff < 0.50) ul = 150.5;
     if(seff >= 0.50 && seff < 0.55) ul = 155.9;
     return ul;
}
float getExpectedP1UpperLimit_T2BW_LM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 161.1;
     if(seff >= 0.05 && seff < 0.10) ul = 161.0;
     if(seff >= 0.10 && seff < 0.15) ul = 166.1;
     if(seff >= 0.15 && seff < 0.20) ul = 173.6;
     if(seff >= 0.20 && seff < 0.25) ul = 181.8;
     if(seff >= 0.25 && seff < 0.30) ul = 191.7;
     if(seff >= 0.30 && seff < 0.35) ul = 203.7;
     if(seff >= 0.35 && seff < 0.40) ul = 216.5;
     if(seff >= 0.40 && seff < 0.45) ul = 228.9;
     if(seff >= 0.45 && seff < 0.50) ul = 242.4;
     if(seff >= 0.50 && seff < 0.55) ul = 257.8;
     return ul;
}
float getExpectedM1UpperLimit_T2BW_LM150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 88.7;
     if(seff >= 0.05 && seff < 0.10) ul = 88.1;
     if(seff >= 0.10 && seff < 0.15) ul = 89.1;
     if(seff >= 0.15 && seff < 0.20) ul = 89.2;
     if(seff >= 0.20 && seff < 0.25) ul = 90.4;
     if(seff >= 0.25 && seff < 0.30) ul = 92.0;
     if(seff >= 0.30 && seff < 0.35) ul = 94.2;
     if(seff >= 0.35 && seff < 0.40) ul = 95.5;
     if(seff >= 0.40 && seff < 0.45) ul = 98.0;
     if(seff >= 0.45 && seff < 0.50) ul = 99.8;
     if(seff >= 0.50 && seff < 0.55) ul = 101.9;
     return ul;
}
