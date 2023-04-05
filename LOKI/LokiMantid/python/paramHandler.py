
import MCPL

class ParamHandler:
  '''Handles parameters (userinput, MCPL metadata, defaults) for mantidpython processing'''

  def __init__(self, mcplFile, userInputs, defaultParams, verbose): #add enforce option to 
    self.mcplFile = mcplFile
    self.mcplMetadata = MCPL.MCPLFile(mcplFile).blobs
    self.mcplMetadata.update((k, v.decode()) for k,v in self.mcplMetadata.items() if isinstance(v, bytes)) #decode MCPL byte strings
    self.params = {key: userInputs[key] if userInputs[key] else self.mcplMetadata.get(key, default)
                   for (key,default) in defaultParams.items()} #(order of precedence is userinput > MCPL metadata > defaults)
    #self.params.update((k, v.decode()) for k,v in self.params.items() if isinstance(v, bytes)) #decode MCPL byte strings

    if verbose:
        self.dumpParams()
        #TODO give warning for parameters differing from the defaults?

  def get(self, key):
    try:
      return float(self.params[key])
    except TypeError:
      return self.params[key]

  def getMcplFile(self):
    return self.mcplFile

  def getMcplMetadata(self):
    return self.mcplMetadata

  def dumpMCPLParams(self):
    print(f'    Parameters read from {self.mcplFile}:')
    for key, value in self.mcplMetadata.items():
      print(f'        {key}: {value}')

  def dumpParams(self):
    print(f'    List of used parameters (order of precedence is userinput > MCPL metadata > defaults):')
    for key, value in self.params.items():
      print(f'        {key}: {value}')