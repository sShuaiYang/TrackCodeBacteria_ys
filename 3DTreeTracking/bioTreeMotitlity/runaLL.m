function [bioTree,allDataOri]=runaLL(bioTree)
bioTree=bioTreeTwoPointTracking(bioTree);
bioTree=velocityTracking(bioTree);
[bioTree,allDataOri]=findLongTrace(bioTree);
end