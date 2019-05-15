from spikeFeatures import listAvailableFeatures, FeatureGroup
from computeSpikeFeatures import pickFeatureGroups

class A(FeatureGroup):
	"""docstring for A"""
	_availableFeatures = ( 'a1', 'a2', 'a3')
	
class B(FeatureGroup):
	_availableFeatures = ( 'b1', 'b2', 'b3')

class C(FeatureGroup):
	_availableFeatures = ( 'c1', 'c2', 'c3')

class C1(FeatureGroup):
	_availableFeatures = ( 'c1', )

class AA(FeatureGroup):
	_availableFeatures = ( 'a1', 'b2', 'b3')



print listAvailableFeatures(__name__)

features = ('a1', 'a2')
print features
print pickFeatureGroups(features, __name__)

features = ('a1', 'b2', 'c1')
print features
l = pickFeatureGroups(features, __name__)
print l 
print [x.__name__ for x in l]

features = ('a1', 'b2', 'c1', 'c2')
print features
l = pickFeatureGroups(features, __name__)
print l 
