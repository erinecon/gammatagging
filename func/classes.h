#include "canvas/Persistency/Common/Assns.h"
#include "lardataobj/RecoBase/Hit.h"
#include "dune/GammaTagging/gammatagging.h"
#include "art/Framework/Principal/Handle.h"
#include <string>

template class art::Assns<string, recob::Hit, void>;
template class art::Assns<recob::Hit, string, void>;
template class art::Wrapper<art::Assns<recob::Hit,string,void>>;
template class art::Wrapper<art::Assns<string, recob::Hit, void>>;
