function val = IntegralFunctionLogP(x,A,zTrue,qTrue,pms,mySuffStats,i,j)

pms(i,j,:) = x;
val=LogJoint(false,A,zTrue,qTrue,pms,mySuffStats);

end
