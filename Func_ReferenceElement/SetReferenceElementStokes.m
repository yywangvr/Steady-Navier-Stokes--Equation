function element = SetReferenceElementStokes(elemV,degreeV,elemP,degreeP)
% element = SetReferenceElementStokes(elemV,degreeV,elemP,degreeP)
%    elemV: type of velocity element (0: quadrilatera, 1: triangles, 11: triangle with bubble function)
%    degreeV: interpolation degree for velocity
%    elemP: type of pressure element
%    degreeP: interpolation degree for pressure


elementV = SetReferenceElement(elemV,degreeV); 

element.elemV = elemV; 
element.degreeV = degreeV; 
element.Xe_ref = elementV.Xe_ref; 
element.nenV = elementV.nen; 
element.ngeom = elementV.nen; 
element.ngaus = elementV.ngaus; 
element.GaussPoints = elementV.GaussPoints; 
element.GaussWeights = elementV.GaussWeights; 
element.N = elementV.N; 
element.Nxi = elementV.Nxi; 
element.Neta = elementV.Neta; 

zgp = element.GaussPoints; 

element.elemP = elemP; 
element.degreeP = degreeP; 
element.NP = ShapeFunc(elemP,degreeP,zgp); 
if elemP == 0
    element.nenP = (degreeP+1)^2; 
else
    element.nenP = (degreeP+1)*(degreeP+2)/2; 
end

