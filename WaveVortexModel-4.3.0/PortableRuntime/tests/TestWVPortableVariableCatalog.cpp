#include "WaveVortexRuntime/generated/WVPortableVariableCatalog.hpp"

#include <cstddef>
#include <cstdlib>
#include <set>
#include <string>

using namespace wavevortex::runtime;

int main() {
  const auto require = [](const bool condition) {
    if (!condition) std::abort();
  };
  static_assert(WVPortableVariableCatalog.size() == 23);
  static_assert(portableVariableCatalogBytes() ==
                sizeof(WVPortableVariableCatalog));
  static_assert(findPortableVariable("u") != nullptr);
  static_assert(findPortableVariable("unknown") == nullptr);
  static_assert(portableVariableMetadata(WVPortableVariable::invalid) ==
                nullptr);

  std::set<std::string> names;
  std::set<std::uint8_t> ordinals;
  for (std::size_t index = 0; index < WVPortableVariableCatalog.size();
       ++index) {
    const auto &variable = WVPortableVariableCatalog[index];
    require(variable.ordinal == index);
    require(static_cast<std::size_t>(variable.identifier) == index);
    require(names.insert(variable.name).second);
    require(ordinals.insert(variable.ordinal).second);
    require(portableVariableMetadata(variable.identifier) == &variable);
    require(findPortableVariable(variable.name) == &variable);
  }

  const auto *Ap = findPortableVariable("Ap");
  require(Ap != nullptr);
  require(Ap->kind == WVPortableVariableKind::coefficient);
  require(Ap->naturalRank == WVPortableNaturalRank::coefficient);
  require(Ap->isComplex);
  require(!Ap->isVariableWithLinearTimeStep);
  require(Ap->samplingMask == portableCoefficientSampling);

  const auto *u = findPortableVariable("u");
  require(u != nullptr);
  require(u->kind == WVPortableVariableKind::field);
  require(u->naturalRank == WVPortableNaturalRank::volume);
  require((u->samplingMask & portableFullGridSampling) != 0);
  require((u->samplingMask & portableFixedVerticalProfileSampling) != 0);
  require((u->samplingMask & portablePositionSampling) != 0);
  require(u->primitiveDependencyMask == 1);
  require(u->movingPrimitiveChannel == 0);
  require(u->netCDFAttributeCount == 1);
  require(std::string(u->netCDFAttribute.name) == "standard_name");

  const auto *rhoBar = findPortableVariable("rho_bar");
  require(rhoBar != nullptr);
  require(rhoBar->naturalRank == WVPortableNaturalRank::vertical);
  require(rhoBar->samplingMask == portableFullGridSampling);
  require(rhoBar->dimensionCount == 1);
  require(std::string(rhoBar->dimensions[0]) == "z");

  const auto *zetaZ = findPortableVariable("zeta_z");
  require(zetaZ != nullptr);
  require(zetaZ->primitiveDependencyMask == (16 | 32));
  require(std::string(zetaZ->netCDFAttribute.name) == "short_name");

  return 0;
}
