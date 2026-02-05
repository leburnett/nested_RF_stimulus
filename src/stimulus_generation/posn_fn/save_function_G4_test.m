function save_function_G4_test(func, param, save_dir, filename)
% SAVE_FUNCTION_G4_TEST  Save position function with fixed ID for testing.
%
%   SAVE_FUNCTION_G4_TEST(FUNC, PARAM, SAVE_DIR, FILENAME) saves a position
%   function to .mat and .pfn files, overwriting any existing function
%   with the same ID. This is a testing variant that always uses a fixed
%   function ID to allow quick iteration.
%
%   INPUTS:
%     func     - 1xM array of frame indices (position function values)
%     param    - Structure containing function parameters:
%                .ID      - Function ID number (fixed in calling function)
%                .frames  - Number of frames in associated pattern
%                .gs_val  - Grayscale value (4 for 4-bit patterns)
%                .type    - 'pfn' (position) or 'afn' (analog)
%                .dur     - Total duration in seconds
%     save_dir - Directory to save the function files
%     filename - Base name for the output files
%
%   OUTPUT FILES:
%     - <ID>_<filename>_G4.mat - MATLAB structure with parameters
%     - fun<ID>.pfn - Binary function file for G4 controller
%
%   DIFFERENCES FROM SAVE_FUNCTION_G4:
%     This test version does not automatically increment the function ID,
%     allowing the same function to be overwritten during iterative testing.
%
%   FILE FORMAT:
%     PFN files use a 512-byte header block followed by 16-bit function data.
%     Frame values are stored as (frame_index - 1) for zero-based indexing.
%
%   See also GENERATE_5FLASH_FUNCTION, SAVE_FUNCTION_G4

arguments
    func (1,:) %{mustBeA(func, ["vector"])} % activate this when >MATLAB2020 is requirement
    param (1,1) %{mustBeA(param, ["struct"])}
    save_dir (1,:)
    filename (1,:)
end

param.func = func; %save full function in param structure

%determine function type
if strcmp(param.type,'pfn')
    func = func-1; % frame array starts at 0
    prefix = 'fun';
elseif strcmp(param.type,'afn')
    assert(max(func) <= 10 && min(func) >= -10, 'input exceeds -10 to 10V range')
    func = ADConvert(func); %convert the analog voltage value to a int16 between (-32768, 32767)
    prefix = 'ao';
else
    error('function type must either be "afn" or "pfn"')
end

%create file name
funcname = [prefix num2str(param.ID, '%04d')];

%save header in the first block
block_size = 512; % all data must be in units of block size
Header_block = zeros(1, block_size);

Header_block(1:4) = dec2char(length(func)*2, 4);     %each function datum is stored in two bytes in the currentFunc card
Header_block(5) = length(funcname);
Header_block(6: 6 + length(funcname) -1) = funcname;

%concatenate the header data with function data
Data = signed_16Bit_to_char(func);     
Data_to_write = [Header_block Data];
param.size = length(Data);

%save .mat file
matFileName = fullfile(save_dir, strcat(num2str(param.ID,'%04d'), '_', filename, '_G4.mat'));

if strcmp(param.type,'pfn')==1
    pfnparam = param;
    save(matFileName, 'pfnparam');
else
    afnparam = param;
    save(matFileName, 'afnparam');
end

%save function file
fid = fopen(fullfile(save_dir, strcat(funcname, '.', param.type)), 'w');
fwrite(fid, Data_to_write(:),'uchar');
fclose(fid);