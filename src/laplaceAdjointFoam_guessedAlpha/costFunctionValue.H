J = 0.5*gSum( volField * ( Foam::pow( y.internalField() - yd.internalField(), 2 ) \
 + lambda * Foam::pow( u.internalField(), 2 ) ) );
