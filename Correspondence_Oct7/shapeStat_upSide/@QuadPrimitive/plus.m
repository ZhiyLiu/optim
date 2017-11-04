function qP = plus(qP1, qP2)

pos = qP1.pos + qP2.pos;
r = qP1.r * qP2.r;
q = QuatProd(qP1.q, qP2.q);
theta = qP1.theta + qP2.theta;
elongation = qP1.elongation * qP2.elongation;

qP = QuadPrimitive(pos,r,elongation,q,theta,false);