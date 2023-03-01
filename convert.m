



BUS=[bus(:,1),bus(:,3),load(:,5),load(:,6),bus(:,4),bus(:,5),bus(:,6),bus(:,8),bus(:,9),bus(:,2),bus(:,7),2*ones(41,1),0.001*ones(41,1)]
GEN=[gen(:,1),gen(:,2),gen(:,3),gen(:,4),gen(:,5),gen(:,6),gen(:,8),gen(:,14),gen(:,16),gen(:,17)]
BRNCH1=[brch(:,1),brch(:,2),brch(:,3),brch(:,4),brch(:,5),brch(:,6),brch(:,7),brch(:,8),zeros(52,1),zeros(52,1),ones(52,1)]
BRNCH2=[tan(:,1),tan(:,2),tan(:,13),tan(:,14),zeros(17,1),tan(:,19),tan(:,20),tan(:,21),tan(:,16),tan(:,17),tan(:,10)]
BRANCH=[BRNCH1;BRNCH2];