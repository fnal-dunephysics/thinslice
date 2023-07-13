#ifndef SLICEPARAMS_H
#define SLICEPARAMS_H

//const int nwires_in_slice = 20;
//const int nslices = 480/nwires_in_slice;

namespace pi{
  const double thinslicewidth = 10; //cm
  const double Eslicewidth = 50; //MeV
  const double plim = 1000;
  const int nthinslices = 20;
  const double Eslicewidth_t = 50; //MeV
  const int nthinslices_t = 20;

  const int reco_nbins = 11; // with an unphysical underflow
  const int true_nbins = 11;
  const int reco_nbins3D = reco_nbins*reco_nbins*reco_nbins;
  const int true_nbins3D = true_nbins*true_nbins*true_nbins;
  //const double reco_bins[reco_nbins+1] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  //const double reco_bins[reco_nbins+1] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  const double reco_bins[reco_nbins+1] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  const double true_bins[true_nbins+1] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  
  //const double reco_KE[reco_nbins] = {1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0};
  const double reco_KE[reco_nbins] = {500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0};
  //const double reco_KE[reco_nbins] = {1000, 950, 850, 800, 750, 700, 650, 600, 550, 500, 450, 350, 0};
  //const double reco_KE[reco_nbins] = {1000, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0};
  //const double reco_KE[reco_nbins] = {1000, 900, 850, 800, 750, 700, 650, 600, 550, 400, 0};
  const double true_KE[true_nbins] = {500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0};
  
  const int nbinse=100;
  const int nbinsthickness = 100;
}

namespace p{
  // useless. delete later
  const double thinslicewidth = 4; //cm
  const int nthinslices = 24;
  const double Eslicewidth = 10; //MeV
  const int nbinse=12;
  //const int nbinsthickness = 100;

  const int reco_nbins = 11; // with an unphysical underflow
  const int true_nbins = 11;
  const int reco_nbins3D = reco_nbins*reco_nbins*reco_nbins;
  const int true_nbins3D = true_nbins*true_nbins*true_nbins;
  const double reco_bins[reco_nbins+1] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  //const double reco_bins[reco_nbins+1] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  const double true_bins[true_nbins+1] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  const double reco_KE[reco_nbins] = {500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0};
  //const double reco_KE[reco_nbins] = {500, 475, 450, 425, 400, 375, 350, 325, 300, 275, 250, 225, 200, 175, 150, 125, 100, 75, 50, 25, 0};
  const double true_KE[reco_nbins] = {500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 0};
}


#endif
