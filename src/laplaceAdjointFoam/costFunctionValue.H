J = 0.5*gSum( volField * ( 
    Foam::pow( y.internalField() - yd.internalField(), 2 ) \
 + lambda * Foam::pow( u.internalField(), 2 ) ) );

// integrate in the boundary patches
forAll(patchesToIntegrate, i)
{
    label patchi = mesh.boundaryMesh().findPatchID(patchesToIntegrate[i]);
    J += -12.0*gSum(
        mesh.magSf().boundaryField()[patchi] * (y.boundaryField()[patchi]));
}
