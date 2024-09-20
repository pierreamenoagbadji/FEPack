function opres = interpolateObj(P, prefix, op_name, dom)
  structLoc = dom.locateInDomain(P);
  elts = dom.elements(structLoc.elements, :);
  coos = structLoc.barycoos;
  opres = cell(size(P, 1), 1);

  for idI = 1:size(P, 1)
    % Initialization
    solcell = load([prefix, '_', num2str(find(dom.IdPoints == elts(idI, 1))), '.mat']);
    op = solcell.(op_name);
    opres{idI} = coos(idI, 1) * op;

    % Update
    for idJ = 2:size(coos, 2)
      solcell = load([prefix, '_', num2str(find(dom.IdPoints == elts(idI, idJ))), '.mat']);
      op = solcell.(op_name);
      opres{idI} = opres{idI} + coos(idI, idJ) * op;
    end
  end

  if (size(P, 1) == 1)
    opres = opres{1};
  end
end