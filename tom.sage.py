

# This file was *autogenerated* from the file tom.sage
from sage.all_cmdline import *   # import sage library

_sage_const_8192 = Integer(8192); _sage_const_3 = Integer(3); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_4 = Integer(4); _sage_const_6 = Integer(6); _sage_const_882513924928040357450309066397220091277646372086261274321936601596635633984923954144833987953525030355789627450964842020045905213190909141068279026196225949585429237261728503592100444669669864353868767722416582928101113854249933605977465180887102288129524044426194007768010389950938670371868270087297861280289222037038196917328635443893936780396624484877644335395624394774212913049813267702784359867981812156080549117151052105851157969153709059893717204375276474518134054712218071996803020719395732929357944271886502281861502225304523408273112840565200387270048446621707956696583329379915354719635474740128176979543416755700538238895767041555995518758577120323308090557361407430248751899094211602218160554079588129235568769551593586101741605079485545639002850701725580259019934589403612001815759960165642899504924749218789029946547807961155389046773646287096703614538921309113172448574050229723530653824449151532460412232936543106594721852464369376951669296531295769246548879893365861040122026786023658939416130455273784000756651152254514757646665490491541287191882325100445899986904549224546256724848052052451550541753893110709589928876178003543346253926954495934574831987166208379756558755763570708848481789965494470453503175671855597849439832807634008599605516630578586291482531874344321677221824689484874059686330283997548008731468506222692361864842911248305046289917123423048074602333646393225367478980565412564784626856775487260776907529621233126252045766341079859096903323353228659101325219859876114168553050747184124904288193056691854176032627599474893068964903735039077297533935708023540090764282159158066625688764617785166508888880708340899995754595706673205056384685712491413378902027052249448903590807253481926968927675928360953258380032644439318371791380735510329958752069930538552742179201108878509297958405733888040838431184059248741251699210740577562261299379268560255672621793557820398126216147311858134145560478809334474922463707295579391127202137978664567824119960049713581518128600861512092547171873019107991934542466618679140720500035258122356252055328712542546818791022341660301941277628021475625473495113328344864628595154636018819502694903034496592550540236812792044977128767254468169275642143204990874557351217257188338336866868777738916609097528737287599930243000690421795901236778410073183640241930191053947341741361129804397928456818064547844776475678110551368984116533623707281070128537161151532765630295213443164897972605211829915013559 = Integer(882513924928040357450309066397220091277646372086261274321936601596635633984923954144833987953525030355789627450964842020045905213190909141068279026196225949585429237261728503592100444669669864353868767722416582928101113854249933605977465180887102288129524044426194007768010389950938670371868270087297861280289222037038196917328635443893936780396624484877644335395624394774212913049813267702784359867981812156080549117151052105851157969153709059893717204375276474518134054712218071996803020719395732929357944271886502281861502225304523408273112840565200387270048446621707956696583329379915354719635474740128176979543416755700538238895767041555995518758577120323308090557361407430248751899094211602218160554079588129235568769551593586101741605079485545639002850701725580259019934589403612001815759960165642899504924749218789029946547807961155389046773646287096703614538921309113172448574050229723530653824449151532460412232936543106594721852464369376951669296531295769246548879893365861040122026786023658939416130455273784000756651152254514757646665490491541287191882325100445899986904549224546256724848052052451550541753893110709589928876178003543346253926954495934574831987166208379756558755763570708848481789965494470453503175671855597849439832807634008599605516630578586291482531874344321677221824689484874059686330283997548008731468506222692361864842911248305046289917123423048074602333646393225367478980565412564784626856775487260776907529621233126252045766341079859096903323353228659101325219859876114168553050747184124904288193056691854176032627599474893068964903735039077297533935708023540090764282159158066625688764617785166508888880708340899995754595706673205056384685712491413378902027052249448903590807253481926968927675928360953258380032644439318371791380735510329958752069930538552742179201108878509297958405733888040838431184059248741251699210740577562261299379268560255672621793557820398126216147311858134145560478809334474922463707295579391127202137978664567824119960049713581518128600861512092547171873019107991934542466618679140720500035258122356252055328712542546818791022341660301941277628021475625473495113328344864628595154636018819502694903034496592550540236812792044977128767254468169275642143204990874557351217257188338336866868777738916609097528737287599930243000690421795901236778410073183640241930191053947341741361129804397928456818064547844776475678110551368984116533623707281070128537161151532765630295213443164897972605211829915013559); _sage_const_825071401847911350104434016317510780278067283201566754357130842744461195005906395199943718051043673223466730721751949974449248930827451854619414374697619260438961946911749222259780365386291751814038486373465100001426528819046274070279571501847524807503058829803292949581020312542244512310826381340725147604205504421489781205421741398803728752958064593794237603834697854344675329642780388241249095996086652843251798333913223849757378836474373313709743648650234857926624244596951216675142492484734159681684027468402207806205439865347969018612132747592285828382232580799515586663970404103616260646140491130336474241553520642020525516089660137732716330093125235058568988907834750945317769702575635573206226618491699434251943292774344986722857013298412092723076929179307523912203128114090534314456523272023692884978468594858168466441251073816641214891282204410271065824121943315600314597749227889573226269102256807459458734770833396972179130110520264774827679636735786580446623637598159822709594268849412087335528825773591001341459368429154598701743523091608321042347244119166311797831282266269359216738748446338689877553411067459022750494469035170847673664139376175016704738727113860401649828339318679003054736839277206870283192227831617539598383681659913506027472064659239184916246405287069652026544516606976680852350358887886206450202813665604368556089115665756191850387879202714813769550164656889996790902837958128457465941454508790612946593532755196596172284379325906458537834759016596544236993315228152015361031327422453823347007209987459925728347755977454611604649086058581101577344837613218863165735163293232182683983328352703931948934018941328642629738690331539423670935877924451188873671016602730573554563515446477428747122842429919747738097859733596295890237973718047460047408090635470616266779559457882525693610773409953808315167993488749629373478589355791418072020121684946874867922541121111773785091209812728948633266300346422009914729246105964646426829796911062733042949990258625136306447831184696881439764159871493270186815550008827682405281004158733425515557787170576476464934246505105516394418797543153382468729141414796208872453491051467183673152632331727482671336127896268657700582305093675403324150839387303020629188315382088652129003636647961202147261937023281222481237476524670557355735034515267765254502192363955049952594249904232853763728215835537753171258068562061271398805505638547424578150686558691852996512332614361427686993333226212556141099 = Integer(825071401847911350104434016317510780278067283201566754357130842744461195005906395199943718051043673223466730721751949974449248930827451854619414374697619260438961946911749222259780365386291751814038486373465100001426528819046274070279571501847524807503058829803292949581020312542244512310826381340725147604205504421489781205421741398803728752958064593794237603834697854344675329642780388241249095996086652843251798333913223849757378836474373313709743648650234857926624244596951216675142492484734159681684027468402207806205439865347969018612132747592285828382232580799515586663970404103616260646140491130336474241553520642020525516089660137732716330093125235058568988907834750945317769702575635573206226618491699434251943292774344986722857013298412092723076929179307523912203128114090534314456523272023692884978468594858168466441251073816641214891282204410271065824121943315600314597749227889573226269102256807459458734770833396972179130110520264774827679636735786580446623637598159822709594268849412087335528825773591001341459368429154598701743523091608321042347244119166311797831282266269359216738748446338689877553411067459022750494469035170847673664139376175016704738727113860401649828339318679003054736839277206870283192227831617539598383681659913506027472064659239184916246405287069652026544516606976680852350358887886206450202813665604368556089115665756191850387879202714813769550164656889996790902837958128457465941454508790612946593532755196596172284379325906458537834759016596544236993315228152015361031327422453823347007209987459925728347755977454611604649086058581101577344837613218863165735163293232182683983328352703931948934018941328642629738690331539423670935877924451188873671016602730573554563515446477428747122842429919747738097859733596295890237973718047460047408090635470616266779559457882525693610773409953808315167993488749629373478589355791418072020121684946874867922541121111773785091209812728948633266300346422009914729246105964646426829796911062733042949990258625136306447831184696881439764159871493270186815550008827682405281004158733425515557787170576476464934246505105516394418797543153382468729141414796208872453491051467183673152632331727482671336127896268657700582305093675403324150839387303020629188315382088652129003636647961202147261937023281222481237476524670557355735034515267765254502192363955049952594249904232853763728215835537753171258068562061271398805505638547424578150686558691852996512332614361427686993333226212556141099)
n = _sage_const_8192 
l = ceil(n / _sage_const_3 )

def toPoly(x):
    ret = []
    
    
    while(x != _sage_const_0 ):

        ret += [x % (_sage_const_1  << l)]
        x = x >> l
    
    return ret

def eval(p, x):
    
    ret = _sage_const_0 
    
    for i in range(len(p)):
        ret += p[i] * (x ** i)
       
    return ret

def point(p):
    
    ret = []
    ret += [p[_sage_const_2 ]]
    x = [_sage_const_2 , -_sage_const_1 , _sage_const_1 , _sage_const_0 ]
    
    for el in x:
        ret += [eval(p, el)]
    return ret

def calcAB(ABpts):
    ret = []
    
    ret += [ABpts[_sage_const_4 ]] #0
    
    ret += [_sage_const_2  * ABpts[_sage_const_0 ] - ABpts[_sage_const_1 ]/_sage_const_6  - ABpts[_sage_const_2 ]/_sage_const_3  + ABpts[_sage_const_3 ] - ABpts[_sage_const_4 ]/_sage_const_2 ]
    
    ret += [- ABpts[_sage_const_0 ] +  ABpts[_sage_const_2 ]/_sage_const_2  +  ABpts[_sage_const_3 ]/_sage_const_2  -  ABpts[_sage_const_4 ]]
    
    ret += [-_sage_const_2  *  ABpts[_sage_const_0 ] + ABpts[_sage_const_1 ]/_sage_const_6  -  ABpts[_sage_const_2 ]/_sage_const_6  - ABpts[_sage_const_3 ]/_sage_const_2  +  ABpts[_sage_const_4 ]/_sage_const_2  ]
    
    ret += [ ABpts[_sage_const_0 ]]
    
    return ret
    

a = _sage_const_882513924928040357450309066397220091277646372086261274321936601596635633984923954144833987953525030355789627450964842020045905213190909141068279026196225949585429237261728503592100444669669864353868767722416582928101113854249933605977465180887102288129524044426194007768010389950938670371868270087297861280289222037038196917328635443893936780396624484877644335395624394774212913049813267702784359867981812156080549117151052105851157969153709059893717204375276474518134054712218071996803020719395732929357944271886502281861502225304523408273112840565200387270048446621707956696583329379915354719635474740128176979543416755700538238895767041555995518758577120323308090557361407430248751899094211602218160554079588129235568769551593586101741605079485545639002850701725580259019934589403612001815759960165642899504924749218789029946547807961155389046773646287096703614538921309113172448574050229723530653824449151532460412232936543106594721852464369376951669296531295769246548879893365861040122026786023658939416130455273784000756651152254514757646665490491541287191882325100445899986904549224546256724848052052451550541753893110709589928876178003543346253926954495934574831987166208379756558755763570708848481789965494470453503175671855597849439832807634008599605516630578586291482531874344321677221824689484874059686330283997548008731468506222692361864842911248305046289917123423048074602333646393225367478980565412564784626856775487260776907529621233126252045766341079859096903323353228659101325219859876114168553050747184124904288193056691854176032627599474893068964903735039077297533935708023540090764282159158066625688764617785166508888880708340899995754595706673205056384685712491413378902027052249448903590807253481926968927675928360953258380032644439318371791380735510329958752069930538552742179201108878509297958405733888040838431184059248741251699210740577562261299379268560255672621793557820398126216147311858134145560478809334474922463707295579391127202137978664567824119960049713581518128600861512092547171873019107991934542466618679140720500035258122356252055328712542546818791022341660301941277628021475625473495113328344864628595154636018819502694903034496592550540236812792044977128767254468169275642143204990874557351217257188338336866868777738916609097528737287599930243000690421795901236778410073183640241930191053947341741361129804397928456818064547844776475678110551368984116533623707281070128537161151532765630295213443164897972605211829915013559 
b = _sage_const_825071401847911350104434016317510780278067283201566754357130842744461195005906395199943718051043673223466730721751949974449248930827451854619414374697619260438961946911749222259780365386291751814038486373465100001426528819046274070279571501847524807503058829803292949581020312542244512310826381340725147604205504421489781205421741398803728752958064593794237603834697854344675329642780388241249095996086652843251798333913223849757378836474373313709743648650234857926624244596951216675142492484734159681684027468402207806205439865347969018612132747592285828382232580799515586663970404103616260646140491130336474241553520642020525516089660137732716330093125235058568988907834750945317769702575635573206226618491699434251943292774344986722857013298412092723076929179307523912203128114090534314456523272023692884978468594858168466441251073816641214891282204410271065824121943315600314597749227889573226269102256807459458734770833396972179130110520264774827679636735786580446623637598159822709594268849412087335528825773591001341459368429154598701743523091608321042347244119166311797831282266269359216738748446338689877553411067459022750494469035170847673664139376175016704738727113860401649828339318679003054736839277206870283192227831617539598383681659913506027472064659239184916246405287069652026544516606976680852350358887886206450202813665604368556089115665756191850387879202714813769550164656889996790902837958128457465941454508790612946593532755196596172284379325906458537834759016596544236993315228152015361031327422453823347007209987459925728347755977454611604649086058581101577344837613218863165735163293232182683983328352703931948934018941328642629738690331539423670935877924451188873671016602730573554563515446477428747122842429919747738097859733596295890237973718047460047408090635470616266779559457882525693610773409953808315167993488749629373478589355791418072020121684946874867922541121111773785091209812728948633266300346422009914729246105964646426829796911062733042949990258625136306447831184696881439764159871493270186815550008827682405281004158733425515557787170576476464934246505105516394418797543153382468729141414796208872453491051467183673152632331727482671336127896268657700582305093675403324150839387303020629188315382088652129003636647961202147261937023281222481237476524670557355735034515267765254502192363955049952594249904232853763728215835537753171258068562061271398805505638547424578150686558691852996512332614361427686993333226212556141099 
ab = a *b

A = toPoly(abs(a))
print(A)
B = toPoly(abs(b))

Apts = point(A)
Bpts = point(B)

print(Apts)

print()

print(Bpts)

#calc a faire avec un autre algo de mult dans l'idée

ABpts = [Apts[i] * Bpts[i] for i in range(len(Apts))]

print(ABpts);
print()
AB = calcAB(ABpts)

print(AB)
print()
ab_tm = eval(AB, _sage_const_2 **l)

print(ab)
print()
print(ab_tm)



























        

