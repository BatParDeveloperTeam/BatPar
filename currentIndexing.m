function indexIII = currentIndexing (IData, Current, Sensitivity)          
% Given the vector of current --> Index for a specified current, covering pulse head and after-pulse relaxation

        indexIII =  ((Current - Sensitivity) <= IData) + (IData <= (Current + Sensitivity)) == 2; % pre-indexing
                indexIII_copy = indexIII;
                for ccc = 1: length(indexIII) - 1
                    if (IData(ccc) == 0) && (indexIII(ccc+1) == 1)
                        indexIII_copy(ccc) = 1;                            % indexing：add the leading edge (pulse head) of the current pulse
                    end
                    if  (indexIII(ccc) == 1) && (IData(ccc+1) == 0)        % Make sure that current drops to ZERO; Only in this case, can the pulse end be considered 
                        for p_index = ccc + 1 : length(indexIII) - 1 
                            if IData(p_index) == 0 && IData(p_index +1) ~= 0
                                    indexIII_copy(ccc + 1 : p_index) = 1;  % indexing：add the following edge of the current pulse, as well as the relaxation after the pulse
                                break
                            end
                        end
                    end 
                end
        indexIII = indexIII_copy;

end