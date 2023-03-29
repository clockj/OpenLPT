#include "ObjectInfo.h"

ObjectInfo::ObjectInfo (const ObjectInfo& object)
{
    _pt_center      = object._pt_center;
    _tri_error      = object._tri_error;
    _match_pos_info = object._match_pos_info;
}

