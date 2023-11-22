var('x,y,z,t,x0,y0,z0,X,Y,Z')
# p = x + x*y*z + 1
# p = x*y - 2*x*y^2 - z + 3*x^2*z - y*z + 2*y^2*z - 2*x*z^2
p = x^2 + y^2 + z^2 - 1
print( p.subs(x=x0+t*X, y=y0+t*Y, z=z0+t*Z).collect(t) )