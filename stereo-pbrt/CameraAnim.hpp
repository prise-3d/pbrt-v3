#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef _PBRT_CAMERA_ANIM_HPP
#define _PBRT_CAMERA_ANIM_HPP

#include "../core/geometry.h"

#include <string>
#include <vector>

namespace pbrt {
/**
 * \file CameraAnim.hpp définitionde la classe d'animation d'une caméra
 * 
 * Format du fichier : 
 * - le fichier doit contenir le nombre de frames à générer, précédé par le mot-clé numberOfFrames ;
 * - chaque position de caméra est indiquée par le mot-clé LookAt,
 *   suivi de la position de la caméra, du point visé et du vecteur up.
 * Remarque 1 : le schéma d'interpolation des positions de caméras utulisé est celui de 
 * CatmullRom. Il implique que (i) le nombre de positions clé soit supérieur ou égale à 4 et
 * (ii) que (nombre de frame - 1) soit divisible par (nombre de positions clé -1).
 * remarque 2 : si le nombre de positions clé est supérieur au nombre de frames à générer, ce
 * dernier est ignoré.
 */
class CameraAnim {
private:
  std::vector <Point3f> pos, look;
  std::vector <Vector3f> up;
  int curvp;// indice du point de vue courant
  int numberOfFrames; // nombre de frames composant la séquence
  int firstFrame, lastFrame;// début et fin de l'intervalle de calcul
public:
  CameraAnim();
  CameraAnim(std::string filename);

  bool load(std::string filename);

  bool getViewpoint(Point3f &at, Point3f &to, Vector3f &up);

  int getCurrentFrameId(){ return curvp; }
  int getNumberOfFrame(){ return pos.size(); }

private:
  void interpolCamerasLocation();
    
};

}

#endif
