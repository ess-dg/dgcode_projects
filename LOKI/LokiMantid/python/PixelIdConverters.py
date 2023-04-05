
class Larmor2022GeometryConverter:
  '''PixelId converter from full loki rear bank geometry to Larmor2022 reduced geometry'''
  numberOfPacks = 28
  tubePerLayer = 2 * numberOfPacks
  tubeLowerCut = 2*2 #Pack 1-2
  tubeUpperCut = 18*2 #Pack 19-28
  reductedTubePerLayer = tubeLowerCut + (tubePerLayer - tubeUpperCut) #number of tubes not present in the reduced geometry (per layer)
  
  def __init__(self, pixelPerStraw):
    self.pixelPerStraw = int(pixelPerStraw)
    self.pixelPerTube = 7 * pixelPerStraw
    self.pixelPerLayer = self.pixelPerTube * self.tubePerLayer

  def getLayer(self, id): #zero indexed
    return id // self.pixelPerLayer

  def getLowerIdCut_newGeom(self, id):
    return self.getLayer(id) * self.pixelPerLayer + self.tubeLowerCut * self.pixelPerTube

  def getUpperIdCut_newGeom(self, id):
     return self.getLayer(id) * self.pixelPerLayer + self.tubeUpperCut * self.pixelPerTube

  def getReductionPixelOffset_newGeom(self, id):
    return self.getLayer(id) * self.reductedTubePerLayer * self.pixelPerTube + self.tubeLowerCut * self.pixelPerTube

  def isPixelInReducedGeom(self, id):
    return id >= self.getLowerIdCut_newGeom(id) and id < self.getUpperIdCut_newGeom(id)


class TubeNumberingConverter:
  '''PixelId converter from geometry with old tube numbering (first layer to last layer, top to bottom) to geom with new tube numbering (top to bottom, layer by layer)'''
  numberOfPacks = 28
  tubePerLayer = 2 * numberOfPacks

  def __init__(self, pixelPerStraw):
    self.pixelPerStraw = pixelPerStraw
    self.pixelPerTube = 7 * pixelPerStraw
    self.pixelPerLayer = self.pixelPerTube * self.tubePerLayer

  def getNewTubeId(self, oldTubeId):
    packId = oldTubeId // 8
    inPackTubeId = oldTubeId % 8 #[0-7]
    layerId = inPackTubeId % 4 #[0-3]
    return (layerId * self.tubePerLayer) + (packId * 2) + inPackTubeId // 4

  def getNewPixelId_zeroBased(self, oldPixelId): #works with zero-based indexing
    localPixelId = oldPixelId % self.pixelPerStraw
    strawId = (oldPixelId % self.pixelPerTube) // self.pixelPerStraw
    oldTubeId = oldPixelId // self.pixelPerTube
    newTubeId = self.getNewTubeId(oldTubeId)
    return (newTubeId * self.pixelPerTube) + (strawId * self.pixelPerStraw) + localPixelId

  def getNewPixelId(self, oldPixelId):
    oldPixelId_zeroBased = oldPixelId - 1
    return self.getNewPixelId_zeroBased(oldPixelId_zeroBased) + 1
