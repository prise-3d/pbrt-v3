#include "CameraAnim.hpp"
#include "Interpolation.hpp"

#include <fstream>
#include <iostream>

namespace pbrt {
  
CameraAnim::CameraAnim(){
  curvp = 0;
}

CameraAnim::CameraAnim(std::string filename){
  curvp = 0;
  numberOfFrames = 0;
  firstFrame = lastFrame = -1;
  if(!load(filename)){
    std::cerr << "CameraAnim : ";
    std::cerr << "erreur chargement fichier " << filename << std::endl;
  }
}

bool CameraAnim::load(std::string filename){
  
  std::ifstream in(filename);
  if(in.is_open()==false) return false;

  std::string s;

  // relecture des différents points de vue présents
  // dans le fichier
  in >> s;
  while(!in.eof()){
    if(s=="numberOfFrames"){
      in >> s;
      numberOfFrames = stoi(s);
      firstFrame = 1; lastFrame = numberOfFrames;
    } else if(s=="LookAt"){
      float x, y, z;
      in >> x >> y >> z;
      pos.push_back(Point3f(x, y, z));
      in >> x >> y >> z;
      look.push_back(Point3f(x, y, z));
      in >> x >> y >> z;
      up.push_back(Vector3f(x, y, z));
    } else if(s=="interval"){
      in >> firstFrame >> lastFrame;
      if((firstFrame < 0) || (firstFrame > lastFrame) || (lastFrame>numberOfFrames)){
	in.close();
	return false;
      }
      curvp = firstFrame-1;
    } else {
      getline(in, s);
    }
    in >> s;
  }
  in.close();

  if(numberOfFrames > pos.size()){// Une interpolation de nouveaux points de vue est nécessaire
    if(pos.size()<4){// nb de pts de vue insuffisant pour le schéma d'interpolation
      std::cout << "erreur : le nb de positions de caméra doit être >=4 pour les interpolations";
      std::cout << std::endl;
      return false;
    }
    if((numberOfFrames-1)%(pos.size()-1)!=0){// interpolation impossible
       std::cout << "erreur : la (duree-1) de la séquence doit être un multiple";
       std::cout << " du (nb de positions - 1)" << std::endl;
       return false;
    }
    interpolCamerasLocation();
  }
  return true;
}

bool CameraAnim::getViewpoint(Point3f &at, Point3f &to, Vector3f &up){
  if(curvp<lastFrame){
    at = pos[curvp];
    to = look[curvp];
    up = this->up[curvp];
    curvp++;
    return true;
  }

  return false;
}


void CameraAnim::interpolCamerasLocation(){

  // Calculer les vecteurs direction de visée qui seront utilisés pour l'interpolation
  std::vector <Vector3f> vis(look.size());
  for(int i=0; i<look.size(); i++){
    vis[i] = look[i] - pos[i];
    vis[i] /= vis[i].Length();
    std::cout << "vis = " << vis[i] << std::endl;
  }
				

  // dupliquer la première et la dernière position de caméra
  // afin d'assurer que les points de vue interpolés passent
  // tous par les positions fournies dans le fichier d'animation
  pos.insert(pos.begin(), pos[0]);
  pos.push_back(pos[pos.size()-1]);
  vis.insert(vis.begin(), vis[0]);
  vis.push_back(vis[look.size()-1]);
  up.insert(up.begin(), up[0]);
  up.push_back(up[up.size()-1]);

  // interpolation des positions
  std::vector <Point3f> tabPos(numberOfFrames);
  std::vector <Vector3f> tabUp(numberOfFrames), tabVis(numberOfFrames);;
  
  CInterpolation::CatmullRom(pos.size(),numberOfFrames,pos,tabPos);
  CInterpolation::CatmullRom(vis.size(),numberOfFrames,vis,tabVis);
  CInterpolation::CatmullRom(up.size(),numberOfFrames,up,tabUp);

  // recopier les positions interpolées dans les tableaux correspondants
  pos = tabPos;
  up = tabUp;
  look.clear();
  for(int i=0; i<tabVis.size(); i++)
    look.push_back(pos[i]+tabVis[i]);
  
  // for(int i=0; i<tabPos.size(); i++){
  //   std::cout << tabVis[i] << std::endl; //<< " " << look[i] << " " << up[i] << std::endl;
  // }
		
}

}// namespace pbrt
