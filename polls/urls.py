from django.conf.urls import url
from rest_framework import routers
from . import views

# router = routers.DefaultRouter()
# router.register('OnlineDock', views.online_docking)

router = routers.DefaultRouter()
router.register('params', views.SubmitParamsViewSet)
router.register('onlinedockparams', views.OnlinedockViewSet) ###
urlpatterns = router.urls

urlpatterns += [
    url(r'^online-dock', views.online_docking),   ### main function
    url(r'^prepare-protein', views.prepare_protein),
    url(r'^first-step', views.first_step),
    url(r'^second-step', views.second_step),
    url(r'^third-step', views.third_step),
    url(r'^fourth-step', views.fourth_step),
    url(r'^fifth-step', views.fifth_step),
]