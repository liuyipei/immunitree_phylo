function hashmap = bioinfohashmap(keys,values)
%BIOINFOHASHMAP creates a hash map
%
%   HASHMAP = BIOINFOHASHMAP(KEYS,VALUES) creates a hash map with KEYS and
%   their associated VALUES.  KEYS is an array of cells containing strings
%   or integers. VALUES is an MxN array of cells, where M is the same
%   length as KEYS, containing strings, integers or a combination of both.
%   HASHMAP is a bioinfoHashMap object which extends the Java HashMap
%   class.
%
%   Example:
%
%   keys = {'K03454' 'M27323' 'M15390' 'AY509259' 'L07625' 'M19499'...
%           'M58410' 'AF075269' 'AF115393' 'AY159322' 'M30931'...
%           'M33262' 'AF103818' 'AY340701' 'AF447763' 'AF334679'};
%   values = {'HIV-1 (Zaire)' [1 2 8]  ;
%             'HIV1-NDK (Zaire)' [1 2 8]  ;
%             'HIV-2 (Senegal)' [1 2 8]  ;
%             'HIV2-MCN13' [1 2 8]  ;
%             'HIV-2UC1 (IvoryCoast)' [1 2 8]  ;
%             'SIVMM251 Macaque' [1 2 8]  ;
%             'SIVAGM677A Green monkey' [1 2 7]  ;
%             'SIVlhoest L''Hoest monkeys' [1 2 7]  ;
%             'SIVcpz Chimpanzees Cameroon' [1 2 8]  ;
%             'SIVmnd5440 Mandrillus sphinx' [1 2 8]  ;
%             'SIVAGM3 Green monkeys' [1 2 7]  ;
%             'SIVMM239 Simian macaque' [1 2 8]  ;
%             'CIVcpzUS Chimpanzee' [1 2 8]  ;
%             'SIVmon Cercopithecus Monkeys' [1 2 8]  ;
%             'SIVcpzTAN1 Chimpanzee' [1 2 8]  ;
%             'SIVsmSL92b Sooty Mangabey' [1 2 8]  ;
%             };
%   map = bioinfohashmap(keys,values);
%   map.get('M58410')
%

% $Revision: 1.1.6.5 $ $Date: 2010/12/22 16:19:15 $
% Copyright 2007-2009 The MathWorks, Inc.

if numel(keys)~= size(values,1)
    error(message('bioinfo:bioinfohashmap:ValuesAndKeysNotEqual'));
end

hashmap = com.mathworks.toolbox.bioinfo.util.BioinfoHashMap(keys,values);
