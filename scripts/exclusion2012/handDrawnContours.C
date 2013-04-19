TGraph* get_T2tt_BDT_expp1(){

   int i = 0;
   float x[5];
   float y[5];

   x[i] =  162.5; y[i++]=0;
   x[i] =  265.0; y[i++]=112.5;
   x[i] =  275.0; y[i++]=137.5;
   x[i] =  250.0; y[i++]=137.5;
   x[i] =  137.5; y[i++]=0;

   TGraph *gr = new TGraph(5,x,y);

   gr->SetLineColor(2);
   gr->SetLineStyle(3);

   return gr;
}


TGraph* get_T2tt_obs(){

   int i = 0;
   float x[8];
   float y[8];

   x[i] =  150.0; y[i++]=0;
   x[i] =  262.5; y[i++]=100.0;
   x[i] =  262.5; y[i++]=105.0;
   x[i] =  255.0; y[i++]=112.5;
   x[i] =  250.0; y[i++]=112.5;
   x[i] =  137.5; y[i++]= 25.0;
   x[i] =  137.5; y[i++]=  0.0;
   x[i] =  150.0; y[i++]=  0.0;

   TGraph *gr = new TGraph(8,x,y);

   gr->SetLineWidth(4);

   return gr;
}

TGraph* get_T2tt_obsp1(){

   int i = 0;
   float x[8];
   float y[8];

   x[i] =  150.0; y[i++]=0;
   x[i] =  262.5+30; y[i++]=100.0+25;
   x[i] =  262.5+30; y[i++]=105.0+25;
   x[i] =  255.0+30; y[i++]=112.5+25;
   x[i] =  250.0+30; y[i++]=112.5+25;
   x[i] =  137.5; y[i++]= 25.0;
   x[i] =  137.5; y[i++]=  0.0;
   x[i] =  150.0; y[i++]=  0.0;

   TGraph *gr = new TGraph(8,x,y);

   gr->SetLineWidth(1);
   gr->SetLineStyle(7);

   return gr;
}

TGraph* get_T2tt_obsm1(){

   int i = 0;
   float x[6];
   float y[6];

   x[i] =  162.5; y[i++]= 12.5;
   x[i] =  262.5; y[i++]=100.0;
   x[i] =  262.5; y[i++]=112.5;
   x[i] =  250.0; y[i++]=112.5;
   x[i] =  157.5; y[i++]= 37.5;
   x[i] =  162.5; y[i++]= 12.5;

   TGraph *gr = new TGraph(6,x,y);

   gr->SetLineWidth(1);
   gr->SetLineStyle(7);

   return gr;
}


TGraph* get_T2tt_exp(){

   int i = 0;
   float x[6];
   float y[6];

   x[i] =  162.5;   y[i++]=12.5;
   x[i] =  250.0; y[i++]=100.0;
   x[i] =  262.5; y[i++]=112.5;
   x[i] =  250.0; y[i++]=112.5;
   x[i] =  162.5; y[i++]= 37.5;
   x[i] =  162.5; y[i++]= 12.5;

   TGraph *gr = new TGraph(6,x,y);

   gr->SetLineWidth(3);
   gr->SetLineColor(2);;

   return gr;
}


TGraph* get_T2tt_expm1(){

   int i = 0;
   float x[8];
   float y[8];

   x[i] =  150.0; y[i++]=0;
   x[i] =  262.5+30; y[i++]=100.0+25;
   x[i] =  262.5+30; y[i++]=105.0+25;
   x[i] =  255.0+30; y[i++]=112.5+25;
   x[i] =  250.0+30; y[i++]=112.5+25;
   x[i] =  137.5; y[i++]= 37.5;
   x[i] =  137.5; y[i++]=  0.0;
   x[i] =  150.0; y[i++]=  0.0;

   TGraph *gr = new TGraph(8,x,y);

   gr->SetLineWidth(1);
   gr->SetLineColor(2);
   gr->SetLineStyle(3);

   return gr;
}
