#ifndef LIBRARY_ACCESS_H
#define LIBRARY_ACCESS_H

#include <string>
#include <vector>

class LibraryAccess{

  public:
    void LoadLibraryFromFile(std::string libraryfile, bool reflected, bool reflT0);
    const std::vector<float>* GetReflT0(size_t Voxel);
    const std::vector<float>* GetReflCounts(size_t Voxel);
    const std::vector<float>* GetCounts(size_t Voxel);
    const std::vector<float>* GetLibraryEntries(int VoxID, bool wantReflected);
    std::vector<int> GetVoxelCoords(int id, double position[3]);
    std::vector< std::vector<double> > PhotonLibraryAnalyzer(double _energy, const int _scint_yield, const double _quantum_efficiency, int _voxel);

    LibraryAccess();

  private:
    std::vector<std::vector<float> > table_;
    std::vector<std::vector<float> > reflected_table_;
    std::vector<std::vector<float> > reflT_table_;

    const double gLowerCorner[3] = {0,0,0};//{2.5, -200, 0};
    const double gUpperCorner[3] = {47,40,90};//{202.5, 200, 500};
    const int gxSteps = 47;//40;
    const int gySteps = 40;//80;
    const int gzSteps = 90;//100;



};



#endif
