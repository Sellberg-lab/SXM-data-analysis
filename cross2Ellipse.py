def cross2Ellipse(crossImage):
    """
    cross2Ellipse takes an image with black crosses and interprets them as the 
    dimention of elliptical particles. A mask image containing the elliptical 
    particles is returned, together with the lengths, width and positions of the 
    ellipses.
    """
    #Binarize image. Set variable lim = 0.1 * 255. Så att alla värden blir svarta.
    #Skapa variabel binImg > lim. Denna listan kommer göra om alla värden större
    #än 25 till svart och mindre till vitt. 
    #Lagrar binImg 
    
    #Hitta alla kors, deras regioner. (Matlab regionprops)
    #Ska ge lista med alla regioner. Döp listan till stats. Regionerna ska 
    #innehålla bild och 
    #boundingbox kordinater. Tex: [[image, boundingbox],[image, boundingbox]]
    stats = []
    stats.append([skimage.measure.regionprops(crossImage)])
    
    
    #Kolla längden på listan för antalet kors
    #Skapa tomma arrays/listor för alla längder, brädder, xPositions, yPositions
    #och vinklar som 
    #sen kommer returnas. Dessa kommer sedan fyllas i forloopen.

    #forloopa alla kros och lägg till datan i listorna. I forloopen så:
        #k, antalet loops, går från 1-> len(stats)
        #image = stats[k].image
        #box = stats[k].box
        #image analyseras med hjälp av funktionen CrossInfo(image). 
        #xPos. Det man får tillbaka från CrossInfo är xPos inom den beskärda 
        #bilden.  
        #För att få xPos för hela bilden måste boundingboxens x-värde läggas till.
        #Detta värdet läggs till i listan för xPositioner. 
        #Samma gäller yPos.
        #Appenda även vinkeln till angles, längden till lengths, bredden till 
        #widths etc så att man för varje värde på k appendar värdena till listorna
        #som returnas.

    return [mask, lengths, widths, xPositions, yPositions, angles]
