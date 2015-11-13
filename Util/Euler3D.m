function [r1, r2, r3] = Euler3D(phi, theta, psi)
% Z Y Z rotation from https://en.wikipedia.org/wiki/Euler_angles
c1 = cos(phi(:))';
c2 = cos(theta(:))';
c3 = cos(psi(:))';

s1 = sin(phi(:))';
s2 = sin(theta(:))';
s3 = sin(psi(:))';

r1 = [
  c3.*c2.*c1 - s3.*s1
  c3.*c2.*s1 + s3.*c1
  -c3.*s2
  ];

r2 = [
  -s3.*c2.*c1-c3.*s1
  -s3.*c2.*s1 + c3.*c1
  s3.*s2
  ];

r3 = [
  c1.*s2
  s1.*s2
  c2
  ];
  
  end