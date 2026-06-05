#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <iostream>
#include <vector>

void showBranch(const char* branchName = "cand_charge",
                const char* filename   = "chained_ParticleTree.root",
                Long64_t maxEvents     = 15) {
    // Open the file
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: could not open file " << filename << "!" << std::endl;
        return;
    }

    // Try to find the tree (loop over keys if necessary)
    TTree *t = nullptr;
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();
        if (obj->InheritsFrom(TTree::Class())) {
            t = (TTree*)obj;
            break;
        }
    }

    if (!t) {
        std::cerr << "Error: no TTree found in " << filename << "!" << std::endl;
        return;
    }

    // Find the branch
    TLeaf *leaf = t->FindLeaf(branchName);
    if (!leaf) {
        std::cerr << "Error: branch " << branchName << " not found!" << std::endl;
        return;
    }

    std::string type = leaf->GetTypeName();
    std::cout << "Branch " << branchName << " has type " << type << std::endl;

    // Scalars
    float   fval;   int ival;   double dval;   char cval;
    // Vectors
    std::vector<float>   *vfloat   = nullptr;
    std::vector<int>     *vint     = nullptr;
    std::vector<double>  *vdouble  = nullptr;
    std::vector<char>    *vchar    = nullptr;

    // Set branch addresses
    if (type == "float")          t->SetBranchAddress(branchName, &fval);
    else if (type == "int")       t->SetBranchAddress(branchName, &ival);
    else if (type == "double")    t->SetBranchAddress(branchName, &dval);
    else if (type == "char")      t->SetBranchAddress(branchName, &cval);
    else if (type == "vector<float>")   t->SetBranchAddress(branchName, &vfloat);
    else if (type == "vector<int>")     t->SetBranchAddress(branchName, &vint);
    else if (type == "vector<double>")  t->SetBranchAddress(branchName, &vdouble);
    else if (type == "vector<char>")    t->SetBranchAddress(branchName, &vchar);
    else {
        std::cerr << "Unsupported branch type: " << type << std::endl;
        return;
    }

    // Loop over entries
    Long64_t nentries = t->GetEntries();
    if (maxEvents < 0 || maxEvents > nentries) maxEvents = nentries;

    for (Long64_t i = 0; i < maxEvents; i++) {
        t->GetEntry(i);
        std::cout << "Event " << i << ": ";

        if (type == "float")  std::cout << fval;
        if (type == "int")    std::cout << ival;
        if (type == "double") std::cout << dval;
        if (type == "char")   std::cout << (int)cval;

        if (vfloat)   for (auto x : *vfloat)   std::cout << x << " ";
        if (vint)     for (auto x : *vint)     std::cout << x << " ";
        if (vdouble)  for (auto x : *vdouble)  std::cout << x << " ";
        if (vchar)    for (auto x : *vchar)    std::cout << (int)x << " ";

        std::cout << std::endl;
    }
}
